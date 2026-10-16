"""Reading and writing the registry file: load-time checks, the lock, atomic saves and transactions.

``prepare.data.registry`` re-exports these names. The functions here look that module up at call
time for the dataclasses, the schema version, the lock settings, the clock and the helper
functions, because it imports this module and so cannot be imported from here at load time.
"""

import json
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional, Union

from prepare.core.constants import SHORT_LOCK_HEARTBEAT_SECONDS, SHORT_LOCK_POLL_SECONDS
from prepare.core.exceptions import DataAccessError
from prepare.core.settings import settings_or
from prepare.data.file_io import write_text_atomic
from prepare.utils.lockfile import LockPolicy, verify_held

if TYPE_CHECKING:
    from prepare.data.registry import Registry

# The top-level keys the Registry dataclass models; any other key is kept in ``Registry.extra``.
KNOWN_KEYS = frozenset({"version", "created", "updated", "genomes", "datasets", "project", "store"})


def parse_registry_document(text: str, target: Path, supported_version: int) -> Dict[str, Any]:
    """The registry JSON ``text`` read from ``target`` as a dict whose ``version`` is an int.

    Raises ``DataAccessError`` for invalid JSON, a document that is not an object, a version that
    is not a number, or a version above ``supported_version`` (a newer prepare wrote the file).
    """
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise DataAccessError(f"Registry is not valid JSON: {target} ({e})") from e
    if not isinstance(data, dict):
        raise DataAccessError(f"Registry is not a JSON object: {target}")
    try:
        version = int(data.get("version", supported_version))
    except (TypeError, ValueError) as e:
        raise DataAccessError(f"Registry version is not a number: {target}") from e
    if version > supported_version:
        raise DataAccessError(
            f"Registry {target} was written by a newer prepare (schema {version}); upgrade prepare or use that version"
        )
    return {**data, "version": version}


def _registry_module() -> ModuleType:
    """The ``prepare.data.registry`` module, looked up at call time because that module imports this one."""
    from prepare.data import registry

    return registry


def load_registry(path: Optional[Union[str, Path]] = None) -> "Registry":
    """Load the registry, or an empty one bound to the path when the file does not exist."""
    r = _registry_module()
    target = r.registry_path(path)
    if not target.exists():
        return r.Registry(path=target)
    try:
        text = target.read_text()
    except OSError as e:
        raise DataAccessError(f"Cannot read registry {target}: {e}") from e
    data = parse_registry_document(text, target, r.SCHEMA_VERSION)
    return r.Registry(
        version=data["version"],
        created=str(data.get("created", r._now())),
        updated=str(data.get("updated", r._now())),
        genomes=dict(data.get("genomes", {})),
        datasets=dict(data.get("datasets", {})),
        project=dict(data.get("project", {})),
        store=dict(data.get("store", {})),
        extra={k: v for k, v in data.items() if k not in KNOWN_KEYS},
        path=target,
    )


@contextmanager
def _acquire_lock(lock: Path) -> Iterator[Path]:
    """Hold the registry lock file ``lock`` for a with-block (``prepare.utils.lockfile``)."""
    r = _registry_module()
    stale, wait = settings_or(registry_lock_stale=r.LOCK_STALE_SECONDS, registry_lock_wait=r.LOCK_WAIT_SECONDS)
    limits = (stale, wait, SHORT_LOCK_POLL_SECONDS, SHORT_LOCK_HEARTBEAT_SECONDS)
    with r.held_lock(lock, LockPolicy("Registry", *limits)) as held:
        yield held


def _write_registry(registry: "Registry", target: Path, lock: Path) -> Path:
    """Write the registry to ``target`` atomically (temp file, compact JSON, rename) while ``lock`` is still held."""
    r = _registry_module()
    verify_held(lock)  # LockLost (a DataAccessError) if a stall let another holder reclaim it
    target.parent.mkdir(parents=True, exist_ok=True)
    registry.updated = r._now()
    registry.path = target
    # An older schema is upgraded on save; fields added by newer versions default to {}.
    registry.version = r.SCHEMA_VERSION
    payload = {
        **registry.extra,
        "version": registry.version,
        "created": registry.created,
        "updated": registry.updated,
        "genomes": registry.genomes,
        "datasets": registry.datasets,
        "project": registry.project,
        "store": registry.store,
    }
    # write_text_atomic removes its temporary file on any exit that did not replace the target,
    # KeyboardInterrupt included, and fsync keeps a power loss from leaving an empty registry.
    try:
        return write_text_atomic(target, json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n", fsync=True)
    except OSError as e:
        raise DataAccessError(f"Cannot write registry {target}: {e}") from e


def save_registry(registry: "Registry", path: Optional[Union[str, Path]] = None) -> Path:
    """Write an already loaded registry atomically under a lock file."""
    r = _registry_module()
    target = Path(path) if path else (registry.path or r.registry_path())
    target.parent.mkdir(parents=True, exist_ok=True)
    with r._acquire_lock(target.with_name(target.name + ".lock")) as lock:
        return r._write_registry(registry, target, lock)


@contextmanager
def registry_transaction(path: Optional[Union[str, Path]] = None) -> Iterator["Registry"]:
    """Load, mutate and save the registry under one lock, so a concurrent edit is never reverted.

    Not re-entrant: nesting a second call to this function (or to `save_registry`) for the
    same registry file inside this block's body raises ``LockReentry`` at once.

    Use this instead of holding one loaded registry across a long run: the file is
    read inside the lock and written back at the end of the block. If the block
    raises, the lock is released and nothing is written.
    """
    r = _registry_module()
    target = r.registry_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with r._acquire_lock(target.with_name(target.name + ".lock")) as lock:
        registry = r.load_registry(target)
        yield registry
        r._write_registry(registry, target, lock)
