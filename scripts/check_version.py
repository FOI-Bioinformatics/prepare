#!/usr/bin/env python3
"""Check that the package version agrees between pyproject.toml and prepare/__init__.py.

With ``--tag vX.Y.Z`` the version must also equal the tag (without the leading ``v``).
Exit status is 1, with a one-line message, on any mismatch; 0 otherwise.
"""

import argparse
import re
import sys
import tomllib
from pathlib import Path
from typing import List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
VERSION_LINE = re.compile(r'^__version__\s*=\s*["\']([^"\']+)["\']', re.MULTILINE)


def read_versions(root: Path = REPO_ROOT) -> tuple:
    """Return (pyproject version, __init__ version) for the checkout at ``root``."""
    with open(root / "pyproject.toml", "rb") as handle:
        project_version = tomllib.load(handle)["project"]["version"]
    match = VERSION_LINE.search((root / "prepare" / "__init__.py").read_text(encoding="utf-8"))
    return project_version, (match.group(1) if match else None)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", help="release tag to compare with, for example v0.10.1")
    args = parser.parse_args(argv)

    project_version, module_version = read_versions()
    if project_version != module_version:
        print(f"ERROR: pyproject.toml version {project_version} differs from prepare.__version__ {module_version}")
        return 1
    if args.tag is not None and args.tag.removeprefix("v") != project_version:
        print(f"ERROR: tag {args.tag} differs from the package version {project_version}")
        return 1
    print(f"Version {project_version} is consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
