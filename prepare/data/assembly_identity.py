"""Identity marker and staged publishing for megahit assembly folders.

An assembly folder ``<ACC>/<genome_id>_assembly`` may carry a marker file, ``.prepare-assembly.json``,
that records the inputs the assembly was built from (the extracted read files by name and size, the
extraction date, the megahit preset, the minimum contig length and any explicit k values). A later
run compares the marker with the inputs it would use now and redoes the assembly only when they
differ. The megahit version is recorded in the marker but is not compared.

megahit writes into a hidden staging folder next to the final one (``unique_temp_path``), and the
result is moved into place by rename (``publish_assembly``), so an interrupted or failed run never
leaves a half-written folder under the final name. Hidden names are ignored by every folder listing
(``visible_files``), so a leftover staging folder is never taken for an assembly; ``sweep_staging``
removes such leftovers and is called while the sample lock is held.
"""

import glob
import json
import logging
import os
import shutil
from dataclasses import dataclass, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from prepare.data.file_io import unique_temp_path, write_text_atomic

logger = logging.getLogger(__name__)

MARKER_NAME = ".prepare-assembly.json"
MARKER_VERSION = 1
CONTIGS_NAME = "final.contigs.fa"
# Middle part of the hidden name an existing folder is moved to while a new one is published; it
# lets ``sweep_staging`` tell that copy from a staging folder and restore it if the publish was cut short.
_ASIDE_TAG = "aside"

AssemblyState = Tuple[str, str]


def _pairs(value: Any) -> Optional[Tuple[Tuple[str, int], ...]]:
    """``[[str, int], ...]`` as a tuple of pairs, or None when the shape does not match."""
    if not isinstance(value, (list, tuple)):
        return None
    pairs: List[Tuple[str, int]] = []
    for item in value:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            return None
        name, number = item
        if not isinstance(name, str) or isinstance(number, bool) or not isinstance(number, int):
            return None
        pairs.append((name, number))
    return tuple(pairs)


def _optional(value: Any, kind: type) -> Tuple[bool, Any]:
    """(valid, value) for a field that is None or of ``kind`` (a bool is not accepted as an int)."""
    if value is None:
        return True, None
    if isinstance(value, bool) or not isinstance(value, kind):
        return False, None
    return True, value


@dataclass(frozen=True)
class AssemblyInputs:
    """The inputs an assembly is built from, as compared between runs.

    ``reads`` holds (file name, size in bytes) pairs of the extracted FASTQ files, sorted by name;
    modification times are not used, since copying a project changes them. ``preset`` is None for
    megahit's own defaults, and ``k_flags`` holds sorted (flag, value) pairs without leading dashes.
    """

    reads: Tuple[Tuple[str, int], ...]
    extraction_date: Optional[str]
    preset: Optional[str]
    min_contig_len: Optional[int]
    k_flags: Tuple[Tuple[str, int], ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        """A JSON-ready dictionary; ``from_dict`` reads it back."""
        return {
            "reads": [[name, size] for name, size in self.reads],
            "extraction_date": self.extraction_date,
            "preset": self.preset,
            "min_contig_len": self.min_contig_len,
            "k_flags": [[flag, value] for flag, value in self.k_flags],
        }

    @classmethod
    def from_dict(cls, data: Any) -> Optional["AssemblyInputs"]:
        """The inputs ``to_dict`` wrote, or None when ``data`` does not have that shape."""
        if not isinstance(data, dict):
            return None
        reads = _pairs(data.get("reads"))
        k_flags = _pairs(data.get("k_flags", []))
        date_ok, extraction_date = _optional(data.get("extraction_date"), str)
        preset_ok, preset = _optional(data.get("preset"), str)
        length_ok, min_contig_len = _optional(data.get("min_contig_len"), int)
        if reads is None or k_flags is None or not (date_ok and preset_ok and length_ok):
            return None
        return cls(reads, extraction_date, preset, min_contig_len, k_flags)

    def differences(self, other: "AssemblyInputs") -> List[str]:
        """One ``field: recorded -> now`` line per field that differs, ``self`` being the recorded side."""
        changed: List[str] = []
        for item in fields(self):
            before, after = getattr(self, item.name), getattr(other, item.name)
            if before != after:
                changed.append(f"{item.name}: {before!r} -> {after!r}")
        return changed


def _file_size(path: Path) -> int:
    """Size of ``path`` in bytes, or -1 when it cannot be read (a missing file then counts as a change)."""
    try:
        return path.stat().st_size
    except OSError:
        return -1


def assembly_inputs(
    reads: Sequence[Union[str, Path]],
    extraction_date: Optional[str],
    preset: Optional[str],
    min_contig_len: Optional[int],
    k_flags: Optional[Mapping[str, int]] = None,
) -> AssemblyInputs:
    """The ``AssemblyInputs`` for a set of extracted read files and megahit settings.

    Read files are identified by name and size, not modification time. A preset of ``"default"``
    becomes None (both leave the flag out), and ``k_flags`` keys lose any leading dashes.
    """
    pairs = sorted((Path(path).name, _file_size(Path(path))) for path in reads)
    normalised_preset = None if preset in (None, "default") else str(preset)
    flags = sorted((str(key).lstrip("-"), int(value)) for key, value in (k_flags or {}).items())
    return AssemblyInputs(tuple(pairs), extraction_date, normalised_preset, min_contig_len, tuple(flags))


def read_marker(out_dir: Union[str, Path]) -> Optional[Dict[str, Any]]:
    """The marker in ``out_dir`` as a dictionary, or None when it is missing or not a JSON object."""
    path = Path(out_dir) / MARKER_NAME
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def write_marker(
    out_dir: Union[str, Path], inputs: AssemblyInputs, version: str, params: Optional[Mapping[str, Any]]
) -> Path:
    """Write the marker for ``inputs`` into ``out_dir`` (atomically) and return its path.

    ``version`` is the megahit version line and ``params`` any further settings worth keeping
    (threads, memory); both are recorded for reference and are not compared by ``assembly_state``.
    """
    marker = {
        "marker_version": MARKER_VERSION,
        "inputs": inputs.to_dict(),
        "megahit_version": version,
        "params": dict(params or {}),
        "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return write_text_atomic(Path(out_dir) / MARKER_NAME, json.dumps(marker, indent=2, sort_keys=True) + "\n")


def assembly_state(
    out_dir: Union[str, Path], expected: Optional[AssemblyInputs], accept_unmarked: bool
) -> AssemblyState:
    """Whether ``out_dir`` holds an assembly of ``expected``, as (state, reason).

    States: ``absent`` (no folder), ``incomplete`` (a folder without ``final.contigs.fa``),
    ``current`` or ``stale``. With ``expected`` None, any folder with contigs is current (the
    inputs are not checked). Otherwise a marker is compared with ``expected``; a folder without a
    marker is current when ``accept_unmarked`` is True and stale otherwise, and an unreadable
    marker is stale.
    """
    path = Path(out_dir)
    if not path.is_dir():
        return "absent", "no assembly folder"
    if not (path / CONTIGS_NAME).is_file():
        return "incomplete", f"no {CONTIGS_NAME}"
    if expected is None:
        return "current", f"{CONTIGS_NAME} present"
    if not (path / MARKER_NAME).exists():
        return ("current", "no marker; accepted as it is") if accept_unmarked else ("stale", "no marker")
    marker = read_marker(path)
    recorded = AssemblyInputs.from_dict(marker.get("inputs")) if marker is not None else None
    if recorded is None:
        return "stale", "marker unreadable"
    changed = recorded.differences(expected)
    if changed:
        return "stale", "inputs changed: " + "; ".join(changed)
    return "current", "marker matches"


def _aside_path(out_dir: Path) -> Path:
    """A hidden name next to ``out_dir`` for the existing folder while a new one is published."""
    return unique_temp_path(out_dir.with_name(f"{out_dir.name}.{_ASIDE_TAG}"))


def publish_assembly(staging: Union[str, Path], out_dir: Union[str, Path]) -> None:
    """Move the finished ``staging`` folder to ``out_dir``, replacing any folder there.

    An existing ``out_dir`` is first renamed to a hidden name beside it; if the rename of
    ``staging`` then fails, that copy is renamed back and the error is raised. If renaming it back
    fails too, that is logged and the first error is raised; the copy stays beside ``out_dir``
    under its hidden name until the caller's cleanup sweeps it, and the next run assembles again
    (the registry record was cleared before the run). Once the new folder is in place the old copy
    is removed.
    """
    source, target = Path(staging), Path(out_dir)
    aside: Optional[Path] = None
    if target.exists() or target.is_symlink():
        aside = _aside_path(target)
        os.rename(target, aside)
    try:
        os.rename(source, target)
    except OSError:
        if aside is not None:
            try:
                os.rename(aside, target)
            except OSError as restore_error:
                logger.error("Could not put %s back as %s: %s", aside.name, target, restore_error)
        raise
    if aside is not None:
        _remove(aside)


def _remove(path: Path) -> None:
    """Remove a folder, a file or a link."""
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path, ignore_errors=True)
    else:
        path.unlink(missing_ok=True)


def _newest(paths: List[Path]) -> Path:
    """The entry of ``paths`` modified last; an entry that cannot be stat'ed counts as oldest."""

    def _mtime(path: Path) -> float:
        try:
            return path.stat().st_mtime
        except OSError:
            return 0.0

    return max(paths, key=_mtime)


def _interrupted_publish_result(asides: List[Path], stagings: List[Path]) -> Path:
    """The folder to put back after a publish stopped between its two renames.

    A staging folder holding both ``final.contigs.fa`` and the marker is the finished new assembly
    (the marker is written last); without one, the old copy moved aside.
    """
    finished = [p for p in stagings if (p / CONTIGS_NAME).is_file() and (p / MARKER_NAME).is_file()]
    return _newest(finished) if finished else _newest(asides)


def sweep_staging(out_dir: Union[str, Path]) -> List[Path]:
    """Remove the hidden ``.<name>.*.tmp`` entries beside ``out_dir`` and return them.

    These are staging folders, or old copies moved aside, that an interrupted run left behind. If
    ``out_dir`` itself is missing because a run stopped between moving the old folder aside and
    renaming the new one in (an old copy and a staging folder are both there), the publish is
    finished first: a staging folder that holds both ``final.contigs.fa`` and the marker is the
    new assembly, written completely, and is renamed into place; otherwise the most recent old
    copy is renamed back. An old copy with no staging folder beside it is what a removal that
    failed after a finished publish leaves; ``out_dir`` missing then means it was removed after
    that publish (by hand, for example), so the old copy is removed rather than brought back. The
    caller holds the sample lock, so no other run is writing to these names.
    """
    target = Path(out_dir)
    if not target.parent.is_dir():
        return []
    leftovers = sorted(Path(p) for p in glob.glob(str(target.parent / f".{glob.escape(target.name)}.*.tmp")))
    asides = [p for p in leftovers if p.name.startswith(f".{target.name}.{_ASIDE_TAG}.")]
    stagings = [p for p in leftovers if p not in asides]
    if asides and stagings and not (target.exists() or target.is_symlink()):
        restored = _interrupted_publish_result(asides, stagings)
        os.rename(restored, target)
        leftovers.remove(restored)
        logger.warning("Restored %s from %s, left by an interrupted run", target, restored.name)
    for path in leftovers:
        _remove(path)
        logger.info("Removed %s, left by an interrupted assembly", path)
    return leftovers
