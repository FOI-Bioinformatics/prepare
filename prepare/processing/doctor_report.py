"""Checks behind ``prepare doctor``: is this environment ready to run PREPARE?

Each check function returns one or more ``Check(name, status, detail, data)`` records, where
``status`` is ``ok``, ``warn`` or ``fail``, ``detail`` is one line for a person and ``data`` holds
the same facts for ``--json``. No check raises for a problem it is looking for: a missing tool,
an unreadable config file, a store without a marker or an unreachable service all become a
``Check``. ``run_checks`` runs them in a fixed order; the command decides the exit code from
``overall_status``.

What counts as a failure:

- a tool older than its floor in ``prepare.utils.tools.TOOLS``, or missing when the command
  named with ``--for`` needs it (missing otherwise is a warning);
- a config file or ``PREPARE_*`` variable that does not parse;
- a store root that is named but has no marker, cannot be reached or is not writable;
- a project registry that cannot be read;
- with ``--network``, NCBI or Branchwater not answering within 10 s.

Free space below ``min_free_gb`` at the project, the temporary folder or the ``.sra`` cache is a
warning, as is an unknown key in the config file.
"""

import logging
import os
import platform
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Union

import requests

from prepare import __version__
from prepare.core import settings
from prepare.core.exceptions import ConfigurationError, DataAccessError
from prepare.data.branchwater_search import DEFAULT_SERVER as BRANCHWATER_SERVER
from prepare.data.registry import load_registry, registry_path
from prepare.data.sra.cleanup import PROJECT_STAGING_FOLDER
from prepare.store.layout import store_paths
from prepare.store.resolve import resolve_store_root
from prepare.utils import resources
from prepare.utils.tools import TOOLS, ToolSpec, ToolStatus, probe_tool

logger = logging.getLogger(__name__)

OK, WARN, FAIL = "ok", "warn", "fail"
_RANK = {OK: 0, WARN: 1, FAIL: 2}

GB = 1024**3
NETWORK_TIMEOUT = 10
NCBI_EINFO_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
SLURM_VARIABLES = (
    "SLURM_JOB_ID",
    "SLURM_ARRAY_TASK_ID",
    "SLURM_CPUS_PER_TASK",
    "SLURM_MEM_PER_NODE",
    "SLURM_MEM_PER_CPU",
    "SLURM_JOB_NODELIST",
)
MIN_PYTHON = (3, 12)


@dataclass(frozen=True)
class Check:
    """The outcome of one check: ``status`` is ``ok``, ``warn`` or ``fail``; ``data`` is for ``--json``."""

    name: str
    status: str
    detail: str
    data: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """The check as a plain dict, as ``doctor --json`` writes it."""
        return {"name": self.name, "status": self.status, "detail": self.detail, "data": dict(self.data)}


def overall_status(checks: Iterable[Check]) -> str:
    """The worst status among ``checks`` (``ok`` when there are none)."""
    worst = OK
    for check in checks:
        if _RANK[check.status] > _RANK[worst]:
            worst = check.status
    return worst


def tools_needed_for(command: str) -> Set[str]:
    """The tools every run of ``command`` needs: those that list it in ``used_by`` and are not optional."""
    return {spec.name for spec in TOOLS.values() if command in spec.used_by and not spec.optional}


# --- Python and tools ------------------------------------------------------------


def check_python() -> Check:
    """The Python interpreter and PREPARE version; fail below Python 3.12."""
    version = platform.python_version()
    data = {"python": version, "executable": sys.executable, "prepare": __version__}
    status = OK if sys.version_info[:2] >= MIN_PYTHON else FAIL
    detail = f"Python {version} ({sys.executable}), PREPARE {__version__}"
    if status == FAIL:
        detail += f"; PREPARE needs Python {'.'.join(map(str, MIN_PYTHON))} or later"
    return Check("python", status, detail, data)


def _tool_data(spec: ToolSpec, status: ToolStatus) -> Dict[str, Any]:
    return {
        "path": status.path,
        "version": status.version_string or None,
        "version_text": status.version_text or None,
        "min_version": spec.min_version,
        "conda_package": spec.conda_package,
        "used_by": list(spec.used_by),
        "optional": spec.optional,
    }


def check_tool(name: str, required: bool = False) -> Check:
    """One tool's path and version against its floor; missing is ``fail`` when ``required``, else ``warn``."""
    spec = TOOLS[name]
    status = probe_tool(name)
    data = _tool_data(spec, status)
    label = f"tool {name}"
    if not status.found:
        use = spec.note or f"used by {', '.join(spec.used_by)}"
        return Check(label, FAIL if required else WARN, f"{status.problem()} ({use})", data)
    if status.not_runnable:
        return Check(label, FAIL if required else WARN, status.problem() or "", data)
    if not status.meets_floor:
        return Check(label, FAIL, status.problem() or "", data)
    floor = f" (needs {spec.min_version} or later)" if spec.min_version else ""
    if status.version is None:
        reason = status.error or "printed no version number"
        return Check(label, WARN, f"{status.path}: version not read ({reason}){floor}", data)
    return Check(label, OK, f"{status.version_text} at {status.path}{floor}", data)


def check_tools(required: Iterable[str] = ()) -> List[Check]:
    """``check_tool`` for every tool in the table, in table order."""
    needed = set(required)
    return [check_tool(name, required=name in needed) for name in TOOLS]


# --- configuration -----------------------------------------------------------------


def _settings_data(resolved: Mapping[str, settings.Resolved]) -> Dict[str, Dict[str, Any]]:
    """``{name: {value, source}}`` for every setting, with secret values hidden."""
    data = {}
    for name, item in resolved.items():
        value = item.value
        if settings.SETTINGS[name].secret and value is not None:
            value = "(set)"
        data[name] = {"value": value, "source": item.source}
    return data


def check_config(args: Any = None, startup_error: Optional[str] = None) -> Check:
    """The config file parses and every setting resolves; the data lists each value and its source.

    ``startup_error`` is the message ``main()`` caught when it could not activate the settings.
    An unknown key in ``[runtime]`` is a warning.
    """
    path = settings.config_path()
    data: Dict[str, Any] = {"path": str(path), "exists": path.exists()}
    if startup_error:
        return Check("config", FAIL, startup_error, data)
    try:
        resolved = settings.resolve_all(args)
        settings.RuntimeSettings.from_resolved(resolved)  # the checks that relate two settings
        warnings = settings.config_warnings(settings.runtime_table())
    except ConfigurationError as e:
        return Check("config", FAIL, str(e), data)
    data["settings"] = _settings_data(resolved)
    changed = [
        f"{name}={item['value']} ({item['source']})"
        for name, item in data["settings"].items()
        if item["source"] != "default"
    ]
    where = str(path) if data["exists"] else f"{path} (not present; built-in defaults)"
    detail = f"{where}; " + (", ".join(changed) if changed else "every setting at its default")
    if warnings:
        data["warnings"] = list(warnings)
        return Check("config", WARN, "; ".join(warnings), data)
    return Check("config", OK, detail, data)


def _setting(name: str) -> Any:
    """The active value of setting ``name``, or its built-in default when the settings do not resolve."""
    try:
        return getattr(settings.active(), name)
    except ConfigurationError:
        return settings.SETTINGS[name].default


# --- store, disk space, registry ------------------------------------------------------


def _existing(path: Path) -> Path:
    """``path`` itself, or its nearest parent that exists (a folder a run would create)."""
    path = Path(path).absolute()
    while not path.exists() and path != path.parent:
        path = path.parent
    return path


def _free_bytes(path: Path) -> Optional[int]:
    """Free bytes on the filesystem holding ``path`` (or its nearest existing parent); None if unreadable."""
    try:
        return int(shutil.disk_usage(_existing(path)).free)
    except OSError as e:
        logger.debug("Could not read free space at %s: %s", path, e)
        return None


def _writable(folder: Path) -> Optional[str]:
    """None when a file can be created and removed in ``folder``, else the error text."""
    try:
        handle, name = tempfile.mkstemp(prefix=".doctor-", dir=folder)
        os.close(handle)
        os.unlink(name)
    except OSError as e:
        return str(e)
    return None


def _registry_store_root(project: Path) -> Optional[str]:
    """The store root the nearest registry records, or None (also when it cannot be read)."""
    try:
        registry = load_registry(registry_path(None, start=project))
    except DataAccessError:
        return None
    root = registry.store.get("root")
    return str(root) if root else None


def check_store(data_root: Optional[str], project: Path, min_free_gb: float) -> Check:
    """The shared data store, if one is configured: marker present, reachable, writable, free space."""
    try:
        root = resolve_store_root(data_root, _registry_store_root(project), require_marker=True)
    except (DataAccessError, ConfigurationError) as e:
        if isinstance(e, ConfigurationError) and str(e).startswith("Config file "):
            # The config check reports this error already; one problem is one failed check.
            return Check("store", WARN, "not checked: the config file does not parse", {"root": data_root})
        return Check("store", FAIL, str(e), {"root": data_root})
    if root is None:
        return Check("store", OK, "no shared data store configured (optional)", {"root": None})
    data: Dict[str, Any] = {"root": str(root), "marker": str(store_paths(root).marker)}
    error = _writable(root)
    data["writable"] = error is None
    if error is not None:
        return Check("store", FAIL, f"store {root} is not writable: {error}", data)
    free = _free_bytes(root)
    data["free_bytes"] = free
    detail = f"store {root}: marker present, writable"
    if free is None:
        return Check("store", WARN, detail + "; free space could not be read", data)
    detail += f", {free / GB:.1f} GB free"
    if min_free_gb and free < min_free_gb * GB:
        return Check("store", WARN, detail + f", below min_free_gb {min_free_gb:g} GB", data)
    return Check("store", OK, detail, data)


def check_free_space(label: str, path: Path, min_free_gb: float) -> Check:
    """Free space where ``path`` would be written; ``warn`` below ``min_free_gb`` (0 turns that off)."""
    name = f"free space {label}"
    free = _free_bytes(path)
    data = {"path": str(path), "free_bytes": free, "min_free_gb": min_free_gb}
    if free is None:
        return Check(name, WARN, f"{path}: free space could not be read", data)
    detail = f"{path}: {free / GB:.1f} GB free"
    if min_free_gb and free < min_free_gb * GB:
        return Check(name, WARN, f"{detail}, below min_free_gb {min_free_gb:g} GB", data)
    return Check(name, OK, detail, data)


def space_locations(project: Path, store_root: Optional[str]) -> List[tuple]:
    """``(label, path)`` for the project, the temporary folder and the ``.sra`` cache.

    The temporary folder is ``temp_folder`` when set, else where a download puts fasterq-dump's
    scratch: the store's ``tmp`` folder, or the project's staging folder ``fastq/.prepare-tmp``.
    """
    store_tmp = store_paths(Path(store_root)).tmp if store_root else None
    temp = _setting("temp_folder") or store_tmp or project / "fastq" / PROJECT_STAGING_FOLDER
    cache = store_tmp or project / "fastq" / ".sra-cache"
    return [("project", project), ("temp", Path(temp)), ("sra-cache", cache)]


def check_registry(project: Path) -> Check:
    """The nearest project registry (walking up from ``project``) loads; none at all is fine."""
    path = registry_path(None, start=project)
    data: Dict[str, Any] = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        return Check("registry", OK, f"no project registry at or above {project} (not a project yet)", data)
    try:
        registry = load_registry(path)
    except DataAccessError as e:
        return Check("registry", FAIL, str(e), data)
    data.update(datasets=len(registry.datasets), genomes=len(registry.genomes))
    detail = f"{path}: {len(registry.datasets)} dataset(s), {len(registry.genomes)} genome(s)"
    return Check("registry", OK, detail, data)


# --- CPUs, memory, scheduler ----------------------------------------------------------


def check_resources(environ: Optional[Mapping[str, str]] = None) -> Check:
    """CPUs available against the node's total, the memory limit found, and any SLURM job variables."""
    env = os.environ if environ is None else environ
    available, total = resources.available_cpus(), os.cpu_count() or 1
    limit = resources.memory_limit_bytes()
    slurm = {name: env[name] for name in SLURM_VARIABLES if env.get(name)}
    data = {"cpus_available": available, "cpus_total": total, "memory_limit_bytes": limit, "slurm": slurm}
    memory = f"memory limit {limit / GB:.1f} GB" if limit else "no memory limit found"
    detail = f"{available} of {total} CPUs available, {memory}"
    if slurm.get("SLURM_JOB_ID"):
        detail += f", SLURM job {slurm['SLURM_JOB_ID']}"
    elif slurm:
        detail += ", SLURM variables set outside a job"
    return Check("resources", OK, detail, data)


# --- network ---------------------------------------------------------------------------


def check_url(name: str, url: str, timeout: float = NETWORK_TIMEOUT) -> Check:
    """One GET to ``url``: any answer below 500 is reachable, 5xx a warning, no answer a failure."""
    data: Dict[str, Any] = {"url": url}
    try:
        response = requests.get(url, timeout=timeout)
    except requests.RequestException as e:
        return Check(name, FAIL, f"{url} not reachable within {timeout:g} s: {e}", data)
    data["status_code"] = response.status_code
    if response.status_code >= 500:
        return Check(name, WARN, f"{url} answered with HTTP {response.status_code}", data)
    return Check(name, OK, f"{url} reachable (HTTP {response.status_code})", data)


def check_network(timeout: float = NETWORK_TIMEOUT) -> List[Check]:
    """NCBI E-utilities (``einfo``) and the Branchwater server, each within ``timeout`` seconds."""
    return [
        check_url("network ncbi", NCBI_EINFO_URL, timeout),
        check_url("network branchwater", BRANCHWATER_SERVER, timeout),
    ]


# --- all of them ------------------------------------------------------------------------


def run_checks(
    project: Union[str, Path] = ".",
    for_command: Optional[str] = None,
    network: bool = False,
    data_root: Optional[str] = None,
    args: Any = None,
    startup_error: Optional[str] = None,
) -> List[Check]:
    """Run every check in order: Python, tools, config, store, free space, registry, resources, network.

    ``for_command`` makes the tools that command needs required (missing is ``fail``); ``args``
    supplies flag values to the settings check; ``startup_error`` is a settings error ``main()``
    already caught.
    """
    project = Path(project).absolute()
    min_free_gb = float(_setting("min_free_gb") or 0)
    checks: List[Check] = [check_python()]
    checks.extend(check_tools(tools_needed_for(for_command) if for_command else ()))
    checks.append(check_config(args, startup_error))
    store = check_store(data_root, project, min_free_gb)
    checks.append(store)
    store_root = store.data.get("root") if store.status != FAIL else None
    checks.extend(check_free_space(label, path, min_free_gb) for label, path in space_locations(project, store_root))
    checks.append(check_registry(project))
    checks.append(check_resources())
    if network:
        checks.extend(check_network())
    return checks


def summary_counts(checks: Sequence[Check]) -> Dict[str, int]:
    """How many checks have each status."""
    counts = {OK: 0, WARN: 0, FAIL: 0}
    for check in checks:
        counts[check.status] += 1
    return counts
