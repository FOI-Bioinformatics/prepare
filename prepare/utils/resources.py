"""CPUs and memory available to this process, as a batch scheduler or container limits them.

``os.cpu_count()`` reports every CPU of the node and ``/proc/meminfo`` all of its memory, while a
SLURM job or a container may use only part of either. These helpers read the limits the kernel
applies (the CPU affinity mask and the cgroup memory limit) and fall back to the SLURM job
variables, so a tool is not started with more threads or memory than the job was given.

The ``/proc`` and ``/sys`` roots are parameters so tests can point them at a fake tree.
"""

import logging
import os
import re
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

# A cgroup v1 limit at or above this is the kernel's way of saying "no limit"
# (9223372036854771712 on 64-bit systems).
_CGROUP_V1_UNLIMITED = 2**60

# Share of a detected memory limit handed to a tool under --assembly-memory auto; the rest is
# left for the Python process itself and the page cache.
AUTO_MEMORY_FRACTION = 0.9

_SIZE_UNITS = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
_SIZE_RE = re.compile(r"^(\d+(?:\.\d+)?|\.\d+)([KMGT]?)$", re.IGNORECASE)


def _positive_int(text: Optional[str]) -> Optional[int]:
    """``text`` as a whole number above zero, or None when it is missing or not one."""
    try:
        value = int(str(text).strip())
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def available_cpus() -> int:
    """CPUs this process may run on.

    The affinity mask (``os.sched_getaffinity``, Linux only) reflects a SLURM cpuset or a
    ``taskset``; without it ``SLURM_CPUS_PER_TASK`` is used, then ``os.cpu_count()``, and 1 when
    nothing is known.
    """
    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is not None:
        try:
            cpus = len(getaffinity(0))
        except OSError:
            cpus = 0
        if cpus > 0:
            return cpus
    slurm = _positive_int(os.environ.get("SLURM_CPUS_PER_TASK"))
    if slurm is not None:
        return slurm
    return os.cpu_count() or 1


def _read_text(path: Path) -> Optional[str]:
    """The stripped content of ``path``, or None when it cannot be read."""
    try:
        return path.read_text().strip()
    except (OSError, UnicodeDecodeError):
        return None


def _cgroup_path(proc_root: Path, controller: Optional[str]) -> Optional[str]:
    """This process's cgroup path from ``/proc/self/cgroup``, or None.

    With ``controller`` None the cgroup v2 line (``0::<path>``) is read; otherwise the cgroup v1
    line whose controller list names it (``4:memory:<path>`` or ``7:cpu,memory:<path>``).
    """
    text = _read_text(proc_root / "self" / "cgroup")
    if text is None:
        return None
    for line in text.splitlines():
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        if (controller is None and parts[0] == "0" and not parts[1]) or controller in parts[1].split(","):
            return parts[2].strip() or "/"
    return None


def _smallest_limit(base: Path, group: str, filename: str) -> Optional[int]:
    """The smallest limit in ``filename`` from ``base``/``group`` up to ``base``; None if none is set.

    A parent group's limit also binds its children, so every level is read, not only the leaf.
    A cgroup v2 ``max`` and a cgroup v1 value at or above 2**60 both mean no limit.
    """
    current = base.joinpath(*[part for part in group.split("/") if part])
    limits = []
    while True:
        value = _positive_int(_read_text(current / filename))
        if value is not None and value < _CGROUP_V1_UNLIMITED:
            limits.append(value)
        if current == base or base not in current.parents:
            break
        current = current.parent
    return min(limits) if limits else None


def memory_limit_bytes(proc_root: Union[str, Path] = "/proc", sys_root: Union[str, Path] = "/sys") -> Optional[int]:
    """Memory this process may use, in bytes, or None when no limit can be found.

    Looked up in order: the cgroup v2 ``memory.max`` and then the cgroup v1
    ``memory.limit_in_bytes``, each the smallest along the path from this process's group (as
    ``/proc/self/cgroup`` names it; a SLURM job on a v1 host sits in ``slurm/uid_N/job_M``) to the
    root; then ``SLURM_MEM_PER_NODE``, and ``SLURM_MEM_PER_CPU`` times the available CPUs (both in
    megabytes). On macOS, and on a Linux host without a limit, the result is None.
    """
    proc_root, sys_root = Path(proc_root), Path(sys_root)
    cgroup = sys_root / "fs" / "cgroup"
    v2_group = _cgroup_path(proc_root, None)
    limit = _smallest_limit(cgroup, v2_group, "memory.max") if v2_group is not None else None
    if limit is None:
        v1_group = _cgroup_path(proc_root, "memory") or "/"
        limit = _smallest_limit(cgroup / "memory", v1_group, "memory.limit_in_bytes")
    if limit is None:
        megabytes = _positive_int(os.environ.get("SLURM_MEM_PER_NODE"))
        if megabytes is None:
            per_cpu = _positive_int(os.environ.get("SLURM_MEM_PER_CPU"))
            megabytes = per_cpu * available_cpus() if per_cpu is not None else None
        limit = megabytes * 1024**2 if megabytes is not None else None
    return limit


def parse_memory(value: str, limit: Optional[int]) -> Optional[Union[int, float]]:
    """Resolve an ``--assembly-memory`` value to what megahit's ``--memory`` receives.

    - ``auto``: ``0.9 * limit`` in bytes when a limit was detected, else None (flag left out).
    - A fraction between 0 and 1 (``0.5``, ``.25``, ``1``): returned as a float; megahit applies
      it to the node's total memory, not to a cgroup limit.
    - A size (``32G``, ``32000M``, ``512K``, ``1T``, binary units) or a whole number of bytes of
      at least 1M: returned as an int.

    Raises:
        ValueError: For anything else, including 0, a number above 1 with a decimal part but no
            unit, and a whole number of bytes below 1M (most likely a size missing its unit).
    """
    text = str(value).strip()
    if text.lower() == "auto":
        return int(AUTO_MEMORY_FRACTION * limit) if limit else None
    match = _SIZE_RE.match(text)
    if match is None:
        raise ValueError(
            f"not a memory value: {value!r} (expected auto, a fraction such as 0.5, or a size such as 32G)"
        )
    number, unit = float(match.group(1)), match.group(2).upper()
    if not unit and 0 < number <= 1:
        return number
    if not unit and "." in match.group(1):
        raise ValueError(f"not a memory value: {value!r} (a fraction is at most 1; give a size with a unit)")
    size = int(number * _SIZE_UNITS[unit])
    if size <= 0:
        raise ValueError(f"not a memory value: {value!r} (must be above zero)")
    if not unit and size < _SIZE_UNITS["M"]:
        raise ValueError(f"not a memory value: {value!r} bytes is too small; give a unit, such as {text}G")
    return size
