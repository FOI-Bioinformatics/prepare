"""
Helpers shared by more than one store command: the timestamp, the no-store hint, the stale
project row, the `.gitignore` guard and the registry's list of linked datasets. The verdict a
store command records for a linked dataset is ``prepare.cli.commands.sra_verdicts.store_verdict``.
"""

import logging
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

from prepare.cli.base import emit_error_json

logger = logging.getLogger(__name__)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _no_store_hint(as_json: bool = False) -> None:
    """Tell the user no store is configured: an ERROR log line, plus a JSON error document on
    stdout with ``--json`` so a script parsing stdout sees it too."""
    message = "No store configured; run: prepare store_init --data-root PATH"
    logger.error(message)
    if as_json:
        emit_error_json(message)


def _stale_project_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """One stale project as reported by ``store_status`` and ``store_gc``.

    Carries the reason staleness was decided ("registry missing" reads very differently from
    "project id differs") and the host that wrote the row, since a project on another
    workstation of a shared store always looks registry-missing from here.
    """
    return {
        "project_id": row["project_id"],
        "name": row.get("name") or row["project_id"],
        "registry": row.get("registry"),
        "hostname": row.get("hostname") or "an unrecorded host",
        "reason": row.get("reason") or "registry missing",
    }


def _gitignore_guard(cwd: Path, log: logging.Logger) -> None:
    """Keep `fastq/` out of git for a project that has just adopted the shared store.

    Only ever reads git state (`git ls-files`) to decide whether to warn; never runs a
    command that changes the git index or working tree. Run from both `store_init` and
    `store_adopt`, since either can be the moment a project's reads become links.
    """
    if not (cwd / ".git").is_dir():
        return

    gitignore = cwd / ".gitignore"
    existing_lines = gitignore.read_text().splitlines() if gitignore.exists() else []
    if not any(line.strip() in ("fastq/", "fastq") for line in existing_lines):
        with gitignore.open("a") as handle:
            if existing_lines and existing_lines[-1] != "":
                handle.write("\n")
            handle.write("fastq/\n")
        log.info("Added fastq/ to %s", gitignore)

    try:
        result = subprocess.run(["git", "ls-files", "fastq"], cwd=cwd, capture_output=True, text=True, check=False)
    except OSError as e:
        log.warning("Could not check git tracking of fastq/: %s", e)
        return

    if result.stdout.strip():
        log.warning("fastq/ is tracked by git; remove it from version control, for example: git rm -r --cached fastq")
