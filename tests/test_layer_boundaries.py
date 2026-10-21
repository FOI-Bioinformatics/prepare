"""
Only the CLI layer imports the CLI layer.

``docs/ARCHITECTURE.md`` keeps ``prepare/cli/`` on top: the data, store, processing, SRA and
plugin layers do not import from it. ``KNOWN_CLI_IMPORTERS`` lists the modules that still do;
it only shrinks.
"""

import ast
from pathlib import Path
from typing import List

PACKAGE_ROOT = Path(__file__).resolve().parent.parent / "prepare"

# The sourmash plugin repeats four CLI commands and borrows two helpers from cli/base.py
# (finding D2 of docs/superpowers/audits/2026-10-21-tech-debt.md).
KNOWN_CLI_IMPORTERS = frozenset({"prepare/plugins/sourmash_plugin.py"})


def _imports_cli(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("prepare.cli"):
            return True
        if isinstance(node, ast.Import) and any(alias.name.startswith("prepare.cli") for alias in node.names):
            return True
    return False


def _cli_importers_outside_the_cli() -> List[str]:
    found = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        relative = path.relative_to(PACKAGE_ROOT.parent).as_posix()
        if relative.startswith("prepare/cli/"):
            continue
        if _imports_cli(path):
            found.append(relative)
    return found


def test_only_the_cli_imports_the_cli():
    unexpected = sorted(set(_cli_importers_outside_the_cli()) - KNOWN_CLI_IMPORTERS)
    assert unexpected == [], f"modules outside prepare/cli/ importing prepare.cli: {unexpected}"


def test_known_cli_importers_still_import_the_cli():
    """An entry that no longer imports the CLI should be removed, so the list keeps shrinking."""
    assert sorted(KNOWN_CLI_IMPORTERS) == _cli_importers_outside_the_cli()
