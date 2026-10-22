"""Size and maintainability ceilings for every module under ``prepare/``.

A module longer than ``MAX_LINES`` lines, or with a radon maintainability index below
``MIN_MAINTAINABILITY``, fails here. Run ``python tests/test_module_sizes.py --check`` to get
the same verdict outside pytest (for example from ``make check``).

``KNOWN_EXCEPTIONS`` lists modules that were already over a ceiling when the guard was added and
are not yet split. Each is held at its size and index of that day, so it can shrink but not grow;
remove its entry once it meets the ceilings.
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

from radon.complexity import cc_visit
from radon.metrics import mi_visit

PACKAGE_ROOT = Path(__file__).resolve().parent.parent / "prepare"
MAX_LINES = 800
MIN_MAINTAINABILITY = 20.0
# Functions at radon rank C (complexity 11-20) sit just under the rank D gate of ``make check``.
# Their number, measured on 2026-10-21, only goes down.
MAX_RANK_C = 104

# path relative to the repository -> (line ceiling, maintainability floor). Empty: every module meets the ceilings.
KNOWN_EXCEPTIONS: Dict[str, Tuple[int, float]] = {}

# Modules split by responsibility in 0.10.1 (Task 6.6); held at this size so they do not regrow.
SPLIT_MODULE_CEILING = 700
SPLIT_MODULES = frozenset(
    {
        "prepare/data/registry.py",
        "prepare/data/registry_io.py",
        "prepare/data/registry_scan.py",
        "prepare/data/read_extraction.py",
        "prepare/data/read_extraction_export.py",
        "prepare/cli/commands/sra.py",
        "prepare/cli/commands/sra_record.py",
        "prepare/data/metadata.py",
        "prepare/data/metadata_fetch.py",
    }
)


def _modules() -> List[Path]:
    return sorted(PACKAGE_ROOT.rglob("*.py"))


def _relative(path: Path) -> str:
    return path.relative_to(PACKAGE_ROOT.parent).as_posix()


def _line_count(path: Path) -> int:
    return len(path.read_text(encoding="utf-8").splitlines())


def _maintainability(path: Path) -> float:
    return float(mi_visit(path.read_text(encoding="utf-8"), multi=True))


def _line_ceiling(relative: str) -> int:
    if relative in SPLIT_MODULES:
        return SPLIT_MODULE_CEILING
    return KNOWN_EXCEPTIONS.get(relative, (MAX_LINES, MIN_MAINTAINABILITY))[0]


def _maintainability_floor(relative: str) -> float:
    return KNOWN_EXCEPTIONS.get(relative, (MAX_LINES, MIN_MAINTAINABILITY))[1]


def oversized_modules() -> List[Tuple[str, int]]:
    """Modules over their line ceiling, as (path relative to the repository, line count)."""
    return [
        (_relative(path), lines) for path in _modules() if (lines := _line_count(path)) > _line_ceiling(_relative(path))
    ]


def unmaintainable_modules() -> List[Tuple[str, float]]:
    """Modules below their maintainability floor, as (path relative to the repository, index)."""
    return [
        (_relative(path), round(score, 1))
        for path in _modules()
        if round(score := _maintainability(path), 1) < _maintainability_floor(_relative(path))
    ]


def test_modules_found():
    assert _modules(), f"no modules found under {PACKAGE_ROOT}"


def test_no_module_exceeds_line_ceiling():
    assert oversized_modules() == []


def test_every_module_meets_maintainability_floor():
    assert unmaintainable_modules() == []


def test_known_exceptions_still_exist_and_are_still_needed():
    """An exception whose module is gone, or already within the ceilings, should be deleted."""
    for relative in KNOWN_EXCEPTIONS:
        path = PACKAGE_ROOT.parent / relative
        assert path.is_file(), f"{relative} no longer exists; remove it from KNOWN_EXCEPTIONS"
        within = _line_count(path) <= MAX_LINES and _maintainability(path) >= MIN_MAINTAINABILITY
        assert not within, f"{relative} now meets the ceilings; remove it from KNOWN_EXCEPTIONS"


def rank_c_blocks() -> List[str]:
    """Every function, method or class at complexity 11-20, as ``path:name``."""
    found = []
    for path in _modules():
        for block in cc_visit(path.read_text(encoding="utf-8")):
            if 11 <= block.complexity <= 20:
                found.append(f"{_relative(path)}:{block.name}")
    return found


def test_rank_c_count_does_not_grow():
    blocks = rank_c_blocks()
    assert len(blocks) <= MAX_RANK_C, f"{len(blocks)} blocks at rank C; the budget is {MAX_RANK_C}:\n" + "\n".join(
        blocks
    )


def test_the_rank_c_budget_is_current():
    """A budget above the count hides progress; lower it when a function leaves rank C."""
    assert len(rank_c_blocks()) == MAX_RANK_C


def test_trivial_code_is_maintainable():
    assert mi_visit("def f():\n    return 1\n", multi=True) >= MIN_MAINTAINABILITY


def main(argv: List[str]) -> int:
    """With ``--check``, report every module over its ceilings and return 1 if there are any."""
    if "--check" not in argv:
        print("usage: python tests/test_module_sizes.py --check")
        return 2
    failures = [f"{path}: {lines} lines (max {_line_ceiling(path)})" for path, lines in oversized_modules()]
    failures += [
        f"{path}: maintainability index {score} (min {_maintainability_floor(path)})"
        for path, score in unmaintainable_modules()
    ]
    for line in failures:
        print(line)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
