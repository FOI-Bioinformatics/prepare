"""
The lint exemption lists in ``setup.cfg`` only shrink.

``per-file-ignores`` lists the package modules that still hold broad ``except Exception``
handlers (``B902``) and those without docstrings on every public name (``D101``-``D103``).
Both lists are meant to be worked down one module at a time; this test freezes their size at
the value measured when it was written, so a new module cannot be added without lowering the
count elsewhere. Lower a constant whenever a module is removed from the list.
"""

import configparser
from pathlib import Path
from typing import List

SETUP_CFG = Path(__file__).resolve().parent.parent / "setup.cfg"

# Measured on 2026-10-21; each only goes down.
B902_EXEMPT = 28
DOCSTRING_EXEMPT = 64


def _package_entries() -> List[str]:
    config = configparser.ConfigParser()
    config.read(SETUP_CFG)
    lines = config["flake8"]["per-file-ignores"].splitlines()
    return [line.strip() for line in lines if line.strip().startswith("prepare/")]


def _count(code: str) -> int:
    return sum(1 for entry in _package_entries() if code in entry.split(":", 1)[1])


def test_broad_except_exemptions_do_not_grow():
    count = _count("B902")
    assert count <= B902_EXEMPT, f"{count} package modules exempt from B902; the list only shrinks ({B902_EXEMPT})"


def test_docstring_exemptions_do_not_grow():
    count = _count("D101")
    assert (
        count <= DOCSTRING_EXEMPT
    ), f"{count} package modules exempt from D101; the list only shrinks ({DOCSTRING_EXEMPT})"


def test_the_constants_are_current():
    """A constant above the measured count hides progress; lower it when a module leaves the list."""
    assert _count("B902") == B902_EXEMPT
    assert _count("D101") == DOCSTRING_EXEMPT
