#!/usr/bin/env python3
"""Check that README.md documents exactly the CLI commands the registry exposes.

Every command name visible in ``prepare --help`` (that is, every registered command whose
``hidden`` property is false; the pre-0.5.0 renamed pointers in ``prepare.cli.commands.renamed``
are ``hidden`` and excluded) must be mentioned somewhere in README.md. Conversely, every literal
``prepare <word>`` invocation shown in README.md must name a command the registry knows about.
A renamed command's old name is allowed to appear as an invocation only inside the "Renamed ...
commands" section that documents the rename; anywhere else it is treated as unknown, the same as
a typo would be.

Usage: python scripts/check_docs_commands.py
Exit status is 1, with the diff printed, when README.md is missing a command or names an unknown
one; 0 when it matches the registry.
"""

import re
import sys
from pathlib import Path
from typing import Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
README_PATH = REPO_ROOT / "README.md"

# Running this file as a script (rather than `python -c` or `-m`) does not add the repository
# root to sys.path, so an unrelated editable install elsewhere (a different checkout or worktree
# pointed at by the active environment's prepare.egg-link) would otherwise shadow the copy of
# prepare this script ships next to. Put the repository root first so the import always resolves
# to this checkout.
sys.path.insert(0, str(REPO_ROOT))

RENAME_HEADING = re.compile(r"^#{2,6}\s*Renamed\b.*\bcommands\b", re.IGNORECASE | re.MULTILINE)
ANY_HEADING = re.compile(r"^#{1,6}\s", re.MULTILINE)
INVOCATION = re.compile(r"^[ \t]*(?:\$ )?prepare\s+([A-Za-z][A-Za-z0-9_-]*)", re.MULTILINE)
# A command table row; the name sits in the first or second cell (the README table has a group column).
TABLE_ROW = re.compile(r"^\|(?:[^|`\n]*\|)?\s*`([A-Za-z][A-Za-z0-9_-]*)`\s*\|", re.MULTILINE)
DOCS_DIR = REPO_ROOT / "docs"


def _registry_commands():
    """Import the CLI package and return its fully-populated command registry."""
    # Imported here, not at module level, so a package import error is reported as this script's
    # own failure rather than as a confusing collection error somewhere else in `make check`.
    from prepare.cli.base import command_registry
    from prepare.cli.main import create_parser

    create_parser()  # populates command_registry via register_all_commands()
    return command_registry.get_all_commands()


def _rename_note_span(text: str) -> Optional[Tuple[int, int]]:
    """Return the (start, end) character offsets of the "Renamed ... commands" section, if any."""
    heading = RENAME_HEADING.search(text)
    if heading is None:
        return None
    following = ANY_HEADING.search(text, heading.end())
    end = following.start() if following else len(text)
    return heading.start(), end


def check(text: str, visible: set, hidden: set) -> Tuple[list, list]:
    """Return (missing command names, unknown (line, word) invocations) for README ``text``."""
    # A command counts as documented only by a table row or an invocation; a passing mention in
    # prose does not.
    documented = {m.group(1) for m in TABLE_ROW.finditer(text)} | {m.group(1) for m in INVOCATION.finditer(text)}
    missing = sorted(name for name in visible if name not in documented)

    rename_span = _rename_note_span(text)
    unknown = []
    for match in INVOCATION.finditer(text):
        word = match.group(1)
        if word in visible:
            continue
        in_rename_note = rename_span is not None and rename_span[0] <= match.start() < rename_span[1]
        if word in hidden and in_rename_note:
            continue
        line_no = text.count("\n", 0, match.start()) + 1
        unknown.append((line_no, word))
    return missing, unknown


def main() -> int:
    if not README_PATH.is_file():
        print(f"ERROR: {README_PATH} not found")
        return 1

    commands = _registry_commands()
    visible = {name for name, command in commands.items() if not command.hidden}
    hidden = {name for name, command in commands.items() if command.hidden}

    text = README_PATH.read_text(encoding="utf-8")
    missing, unknown = check(text, visible, hidden)

    doc_unknown = []
    for page in sorted(DOCS_DIR.glob("*.md")):
        page_text = page.read_text(encoding="utf-8")
        for line_no, word in check(page_text, visible, hidden)[1]:
            doc_unknown.append((page.relative_to(REPO_ROOT).as_posix(), line_no, word))

    if not missing and not unknown and not doc_unknown:
        print(f"README.md documents all {len(visible)} commands, no unknown command named in README.md or docs/")
        return 0

    if missing:
        print(f"README.md is missing {len(missing)} command(s) (a table row or an invocation is required):")
        for name in missing:
            print(f"  - {name}")
    if unknown:
        print(f"README.md names {len(unknown)} command(s) the registry does not know about:")
        for line_no, word in unknown:
            print(f"  README.md:{line_no}: prepare {word}")
    if doc_unknown:
        print(f"docs/ names {len(doc_unknown)} command(s) the registry does not know about:")
        for name, line_no, word in doc_unknown:
            print(f"  {name}:{line_no}: prepare {word}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
