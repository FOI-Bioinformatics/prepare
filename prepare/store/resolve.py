"""
Store root discovery and the small user config file that can name one.

The shared data store is one folder holding a single copy of every
downloaded metagenome, used across organism projects. Its root is found, in
order, from: an explicit flag, the PREPARE_DATA environment variable, a
project registry's recorded root, or the user config file. This module
resolves that precedence and reads/writes the user config.
"""

import json
import logging
import os
from pathlib import Path
from typing import Optional

from prepare.core.constants import STORE_ENV
from prepare.core.exceptions import ConfigurationError, DataAccessError
from prepare.core.settings import config_path, read_config
from prepare.data.file_io import write_text_atomic
from prepare.store.layout import StorePaths, read_marker, store_paths
from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)

# The user config file is shared with the runtime settings, which own its reader; both names
# stay importable from here for existing callers.
__all__ = [
    "config_path",
    "read_config",
    "resolve_optional_store",
    "resolve_store_root",
    "write_config_data_root",
]


def _strip_store_table(text: str) -> str:
    """Remove an existing top-level [store] table block from config text."""
    lines = text.splitlines(keepends=True)
    kept = []
    in_store_table = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("["):
            in_store_table = stripped == "[store]"
            if in_store_table:
                continue
        if in_store_table:
            continue
        kept.append(line)
    return "".join(kept)


def write_config_data_root(root: Path) -> Path:
    """Write the store root into the user config's [store] table.

    Any other top-level tables already present in the file are preserved
    textually; an existing [store] table is replaced.
    """
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)

    existing_text = path.read_text() if path.exists() else ""
    remaining_text = _strip_store_table(existing_text)

    if remaining_text and not remaining_text.endswith("\n"):
        remaining_text += "\n"

    root_posix = Path(root).as_posix()
    # json.dumps produces a double-quoted string with '"' and '\' escaped the
    # same way TOML basic strings require, so it doubles as a safe TOML
    # string literal for a plain path (no control characters to worry about).
    new_block = f"[store]\ndata_root = {json.dumps(root_posix)}\n"

    write_text_atomic(path, remaining_text + new_block, fsync=True)
    return path


class _ConfigRoot:
    """Stands for the config file's ``store.data_root`` until the earlier rules have missed."""

    @staticmethod
    def value() -> Optional[str]:
        store = read_config().get("store", {})
        if not isinstance(store, dict):
            raise ConfigurationError(f"Config file {config_path()}: [store] must be a table")
        root = store.get("data_root")
        if root is not None and not isinstance(root, str):
            raise ConfigurationError(f"Config file {config_path()}: store.data_root must be a string")
        return root


def resolve_store_root(
    explicit: Optional[str],
    registry_root: Optional[str],
    require_marker: bool = True,
) -> Optional[Path]:
    """Resolve the store root, in order: explicit, PREPARE_DATA, registry_root, config.

    Returns None when none of the rules name a root. When a rule names a root
    whose store marker is missing and require_marker is True, raises
    DataAccessError naming the rule that produced the root. On success, the
    resolved root is registered with SecureSubprocess.add_allowed_root.
    """
    candidates = [
        ("--data-root", explicit),
        (STORE_ENV, os.environ.get(STORE_ENV)),
        ("registry store.root", registry_root),
        ("config data_root", _ConfigRoot()),
    ]

    for rule, candidate in candidates:
        if isinstance(candidate, _ConfigRoot):
            candidate = candidate.value()  # the config file is read only when nothing earlier matched
        if not candidate:
            continue

        root = Path(candidate).expanduser().resolve()

        if require_marker and read_marker(root) is None:
            raise DataAccessError(
                f"Store root '{root}' from {rule} has no store marker " "(missing or invalid prepare_store.json)"
            )

        logger.info("Resolved store root from %s: %s", rule, root)
        SecureSubprocess.add_allowed_root(root)
        return root

    return None


def resolve_optional_store(explicit: Optional[str], registry_root: Optional[str]) -> Optional["StorePaths"]:
    """The store's layout for a step that can run without one, or None.

    Reporting and analysis steps read the store when it is there and work from the project's
    own folders when it is not, so an unreachable one (an unmounted volume, a root that has
    moved, a marker that cannot be read) must not stop them: the reason is logged and the step
    continues storeless. Commands that exist to operate on the store itself, and the download
    path that writes into it, keep letting the error through.
    """
    try:
        root = resolve_store_root(explicit, registry_root)
    except DataAccessError as e:
        logger.warning("store unavailable: %s; continuing without it", e)
        return None
    if root is None:
        return None
    return store_paths(root)
