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

from radon.metrics import mi_visit

PACKAGE_ROOT = Path(__file__).resolve().parent.parent / "prepare"
MAX_LINES = 800
MIN_MAINTAINABILITY = 20.0

# path relative to the repository -> (line ceiling, maintainability floor), measured 2026-09-25.
# Follow-up after 0.5.0: split each of these three into a package, as prepare/data/sra.py was
# split in 0.5.0, and drop its entry here once every resulting module meets the ceilings.
# read_extraction.py was re-measured on 2026-09-30, after its assembly code moved to data/assembly.py.
KNOWN_EXCEPTIONS: Dict[str, Tuple[int, float]] = {
    "prepare/data/read_extraction.py": (964, 29.6),
    "prepare/data/registry.py": (1001, MIN_MAINTAINABILITY),
}


def _modules() -> List[Path]:
    return sorted(PACKAGE_ROOT.rglob("*.py"))


def _relative(path: Path) -> str:
    return path.relative_to(PACKAGE_ROOT.parent).as_posix()


def _line_count(path: Path) -> int:
    return len(path.read_text(encoding="utf-8").splitlines())


def _maintainability(path: Path) -> float:
    return float(mi_visit(path.read_text(encoding="utf-8"), multi=True))


def _line_ceiling(relative: str) -> int:
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
