"""``environment.yml``'s version floors match ``prepare.utils.tools.TOOLS`` exactly.

Parsed with a regex rather than a YAML library, since PREPARE has no YAML dependency and this
file's shape (a flat top-level ``dependencies:`` list, one package per line) does not need one.
``TOOLS`` is the single source of truth for which tool needs which floor (``prepare doctor`` and
every command's ``require_tools`` check read it); this test is the guard against the two files
drifting apart silently.
"""

import re
from pathlib import Path
from typing import Dict, Optional

from prepare.utils.tools import TOOLS

REPO_ROOT = Path(__file__).resolve().parent.parent
ENVIRONMENT_YML = REPO_ROOT / "environment.yml"

# A top-level dependency item: exactly two spaces of indent then "- name" or "- name>=1.2.3",
# with an optional trailing "# comment" stripped separately. Excludes the nested "pip:" list
# (indented further) and the "pip:" block header itself.
_DEP_LINE = re.compile(r"^  - ([A-Za-z0-9][A-Za-z0-9_.-]*)(?:(>=)([0-9]+(?:\.[0-9]+)*))?\s*$")


def _parsed_dependencies() -> Dict[str, Optional[str]]:
    """``{package_name: floor_or_None}`` for every top-level conda dependency in environment.yml.

    Only lines under the ``dependencies:`` key are considered, so the ``channels:`` list above it
    (also a two-space-indented ``- name`` list) is not mistaken for a dependency.
    """
    deps: Dict[str, Optional[str]] = {}
    in_dependencies = False
    for raw_line in ENVIRONMENT_YML.read_text(encoding="utf-8").splitlines():
        if raw_line.strip() == "dependencies:":
            in_dependencies = True
            continue
        if not in_dependencies:
            continue
        line = raw_line.split("#", 1)[0].rstrip()
        match = _DEP_LINE.match(line)
        if match is None:
            continue
        name, _op, floor = match.groups()
        deps[name] = floor
    return deps


def _tool_floors_by_conda_package() -> Dict[str, Optional[str]]:
    """``{conda_package: min_version}`` from ``TOOLS``; every tool sharing a package agrees."""
    floors: Dict[str, Optional[str]] = {}
    for spec in TOOLS.values():
        existing = floors.get(spec.conda_package, "unset")
        assert existing in (
            "unset",
            spec.min_version,
        ), f"TOOLS entries sharing conda package {spec.conda_package!r} disagree on min_version"
        floors[spec.conda_package] = spec.min_version
    return floors


def test_environment_yml_parses_at_least_the_known_tool_packages():
    """A sanity check on the regex itself: every installed TOOLS conda package is a dependency line.

    seqkit is the one tool deliberately left commented out (see
    ``test_seqkit_stays_commented_out_and_unfloored``), so it is excluded here.
    """
    deps = _parsed_dependencies()
    for spec in TOOLS.values():
        if spec.name == "seqkit":
            continue
        assert spec.conda_package in deps, f"{spec.conda_package} (provides {spec.name}) missing from environment.yml"


def test_pinned_floors_match_tools_table():
    """Every environment.yml floor for a known tool package equals ``TOOLS[name].min_version``."""
    deps = _parsed_dependencies()
    floors = _tool_floors_by_conda_package()
    mismatches = []
    for package, expected_floor in floors.items():
        actual_floor = deps.get(package)
        if actual_floor != expected_floor:
            mismatches.append(f"{package}: environment.yml has {actual_floor!r}, TOOLS wants {expected_floor!r}")
    assert not mismatches, "\n".join(mismatches)


def test_pigz_and_datasets_floors_are_pinned_explicitly():
    """pigz and ncbi-datasets-cli both carry an explicit floor now, in TOOLS and environment.yml.

    Both used to have ``min_version=None`` in ``TOOLS`` (any version accepted) while
    ``environment.yml`` already pinned ``pigz>=2.4``/``ncbi-datasets-cli>=16`` as documentation
    only: neither floor was actually enforced by ``require_tools``/``prepare doctor``. Named
    explicitly here (on top of the generic comparison in ``test_pinned_floors_match_tools_table``)
    so a future change that quietly drops either floor back to ``None`` fails immediately with a
    clear message, rather than only via the generic mismatch loop above.
    """
    assert TOOLS["pigz"].min_version == "2.4"
    assert TOOLS["datasets"].min_version == "16"
    deps = _parsed_dependencies()
    assert deps["pigz"] == "2.4"
    assert deps["ncbi-datasets-cli"] == "16"


def test_seqkit_stays_commented_out_and_unfloored():
    """seqkit is optional and has no TOOLS floor; it is documented, not installed by default."""
    assert TOOLS["seqkit"].min_version is None
    assert TOOLS["seqkit"].optional is True
    text = ENVIRONMENT_YML.read_text(encoding="utf-8")
    assert re.search(r"^\s*#\s*-\s*seqkit\b", text, re.MULTILINE) is not None
    assert "seqkit" not in _parsed_dependencies()
