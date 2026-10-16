"""Tests for scripts/check_docs_commands.py."""

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "check_docs_commands.py"
spec = importlib.util.spec_from_file_location("check_docs_commands", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

VISIBLE = {"alpha", "beta", "gamma"}


def test_prose_mention_does_not_document_a_command():
    missing, _ = module.check("The alpha step is described here.\n", VISIBLE, set())
    assert "alpha" in missing


def test_table_row_documents_a_command():
    missing, _ = module.check("| `alpha` | does a thing |\n", VISIBLE, set())
    assert "alpha" not in missing


def test_invocations_document_commands_even_indented_or_after_prompt():
    text = "prepare alpha --x\n    prepare beta\n$ prepare gamma\n"
    missing, unknown = module.check(text, VISIBLE, set())
    assert missing == []
    assert unknown == []


def test_unknown_invocation_is_reported_with_line():
    _, unknown = module.check("text\n  $ prepare bogus\n", VISIBLE, set())
    assert unknown == [(2, "bogus")]


def test_hidden_name_allowed_only_in_rename_section():
    text = "## Renamed commands\n\nprepare old\n\n## Other\n\nprepare old\n"
    _, unknown = module.check(text, VISIBLE, {"old"})
    assert unknown == [(7, "old")]


def test_table_row_with_group_column_documents_a_command():
    missing, _ = module.check("| Group | `beta` | purpose |\n| | `gamma` | purpose |\n", VISIBLE, set())
    assert "beta" not in missing and "gamma" not in missing
