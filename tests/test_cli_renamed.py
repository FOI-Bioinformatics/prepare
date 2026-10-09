"""The four SRA analysis commands merged into sra_profile and sra_report in 0.5.0.

Their old names still parse, so a script that calls one gets a pointer to the new command
instead of argparse's list of valid choices, and exit status 2.
"""

import pytest

from prepare.cli.main import create_parser, main

RENAMED = {
    "sra_stats": "sra_profile",
    "sra_profile_quality": "sra_profile",
    "sra-profile-quality": "sra_profile",
    "sra_dashboard": "sra_report",
    "sra-dashboard": "sra_report",
    "sra_compare": "sra_report",
    "sra-compare": "sra_report",
}


@pytest.mark.parametrize("old,new", sorted(RENAMED.items()))
def test_old_name_exits_2_with_a_pointer_on_stderr(old, new, capsys):
    assert main([old]) == 2
    captured = capsys.readouterr()
    assert f"ERROR - prepare.cli.commands.renamed - {old} was renamed: use prepare {new}" in captured.err
    assert captured.out == ""


def test_old_name_with_its_old_flags_still_gets_the_pointer(capsys):
    assert main(["sra_stats", "--fastq-folder", "fastq", "--accessions", "SRR1"]) == 2
    assert "sra_stats was renamed: use prepare sra_profile" in capsys.readouterr().err


def test_registered_renamed_commands_match_the_table():
    from prepare.cli.commands.renamed import RENAMED_COMMANDS

    assert RENAMED_COMMANDS == RENAMED


def test_old_names_are_not_listed_in_the_main_help():
    help_text = create_parser().format_help()
    for old in RENAMED:
        assert old not in help_text
    assert "sra_profile" in help_text and "sra_report" in help_text
