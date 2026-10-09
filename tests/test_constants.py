"""``prepare/core/constants.py``: nothing in the tree still names the constants removed from it.

``DEFAULT_MEMORY_LIMIT_GB``, ``MAX_FILE_SIZE_MB``, ``DEFAULT_PLUGIN_TIMEOUT``, ``ERROR_MESSAGES``
and ``SUCCESS_MESSAGES`` were never read anywhere outside their own definitions, and
``DEFAULT_MAX_WORKERS`` stopped being read when the worker count moved to the CPU detection; this guards
against one being reintroduced (or half-removed, with a stray reference left behind) later.
"""

from pathlib import Path

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parent.parent

REMOVED_CONSTANTS = (
    "DEFAULT_MEMORY_LIMIT_GB",
    "MAX_FILE_SIZE_MB",
    "DEFAULT_PLUGIN_TIMEOUT",
    "ERROR_MESSAGES",
    "SUCCESS_MESSAGES",
    "DEFAULT_MAX_WORKERS",
)

# Source trees that ship with the package or its tests; docs and the plan/report archive under
# docs/superpowers are history, not code, and are allowed to keep mentioning a removed name.
SEARCHED_DIRS = ("prepare", "tests", "scripts")


def test_removed_constants_are_not_referenced_anywhere():
    hits = []
    for dirname in SEARCHED_DIRS:
        for path in (REPO_ROOT / dirname).rglob("*.py"):
            if path == THIS_FILE:
                continue
            text = path.read_text(encoding="utf-8")
            for name in REMOVED_CONSTANTS:
                if name in text:
                    hits.append(f"{path.relative_to(REPO_ROOT)}: {name}")
    assert not hits, "removed constant(s) still referenced:\n" + "\n".join(hits)


def test_removed_constants_are_gone_from_the_module():
    text = (REPO_ROOT / "prepare" / "core" / "constants.py").read_text(encoding="utf-8")
    for name in REMOVED_CONSTANTS:
        assert name not in text
