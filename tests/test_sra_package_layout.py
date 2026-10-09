"""Guards for the layout of the ``prepare.data.sra`` package.

The package re-exports the public names other packages import, and every string patch target
in the tests must name an attribute that exists on the module it names, so a patch never
silently targets a name the code no longer looks up there.
"""

import pathlib
import pkgutil
import re
from typing import List

import prepare.data.sra as sra_package

PATCH_TARGET = re.compile(r'(?:patch|setattr)\(\s*"(prepare\.data\.sra(?:\.[A-Za-z_]+)+)"')


def _bad_patch_targets(text: str) -> List[str]:
    """Patch targets in ``text`` that do not resolve, or that patch a re-export on the package.

    A package-level patch (``prepare.data.sra.<name>``) replaces only the re-exported name; the
    submodules call each other through the owning module, so no internal caller would see it.
    """
    problems = []
    for target in PATCH_TARGET.findall(text):
        owner_name, _, _ = target.rpartition(".")
        try:
            pkgutil.resolve_name(target)
            owner = pkgutil.resolve_name(owner_name)
        except (ImportError, AttributeError) as e:
            problems.append(f"{target}: does not resolve ({e})")
            continue
        if owner is sra_package:
            problems.append(f"{target}: patches the package re-export, not the submodule that owns it")
    return problems


def test_public_names_are_importable_from_the_package():
    from prepare.data import sra

    for name in (
        "download_sra",
        "download_accession",
        "fastq_files",
        "count_fastq_reads",
        "iter_fastq_records",
        "verify_download",
        "parse_verdict_message",
        "compress_fastq",
        "transient_bytes",
        "is_transient_folder",
        "default_max_workers",
        "STORE_READY_STATES",
        "STORE_LINKED_PREFIX",
        "MATE1_SUFFIXES",
        "STOP",
    ):
        assert hasattr(sra, name), name


def test_patch_targets_in_tests_resolve_on_the_owning_submodule():
    problems = []
    for test_file in sorted(pathlib.Path("tests").glob("*.py")):
        if test_file.name == pathlib.Path(__file__).name:
            continue  # this file holds deliberately bad targets for the test below
        problems.extend(f"{test_file}: {p}" for p in _bad_patch_targets(test_file.read_text()))
    assert problems == []


def test_patch_target_guard_catches_setattr_deep_and_package_level_targets():
    text = "\n".join(
        [
            'monkeypatch.setattr("prepare.data.sra.download.no_such_name", fake)',
            'patch("prepare.data.sra.accession.subprocess.no_such_call")',
            'patch(\n    "prepare.data.sra.download_sra")',
            'patch("prepare.data.sra.download.download_sra")',
        ]
    )
    problems = _bad_patch_targets(text)
    assert len(problems) == 3
    assert "no_such_name" in problems[0]
    assert "no_such_call" in problems[1]
    assert "prepare.data.sra.download_sra: patches the package re-export" in problems[2]
