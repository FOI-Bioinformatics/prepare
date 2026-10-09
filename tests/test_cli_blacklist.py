"""Tests for the blacklist command."""

import argparse
import json

from prepare.cli.commands.blacklist import BlacklistCommand


def _args(tmp_path, **kw):
    base = dict(
        add=None,
        remove=None,
        from_file=None,
        reason=None,
        list=False,
        blacklist_file=str(tmp_path / "blacklist.txt"),
        registry=str(tmp_path / "prepare_registry.json"),
    )
    base.update(kw)
    return argparse.Namespace(**base)


def test_add_records_reason_and_writes_file(tmp_path, capsys):
    rc = BlacklistCommand().execute(_args(tmp_path, add=["SRR2517418"], reason="16S amplicon mislabelled as WGS"))
    assert rc == 0
    data = json.loads((tmp_path / "prepare_registry.json").read_text())
    exclusion = data["datasets"]["SRR2517418"]["exclusion"]
    assert exclusion["excluded"] is True and exclusion["reason"] == "16S amplicon mislabelled as WGS"
    assert (tmp_path / "blacklist.txt").read_text().splitlines() == ["SRR2517418  # 16S amplicon mislabelled as WGS"]


def test_add_requires_reason(tmp_path):
    assert BlacklistCommand().execute(_args(tmp_path, add=["SRR1"])) == 1


def test_remove_clears_and_rewrites_file(tmp_path):
    BlacklistCommand().execute(_args(tmp_path, add=["SRR1", "SRR2"], reason="test"))
    assert BlacklistCommand().execute(_args(tmp_path, remove=["SRR1"])) == 0
    data = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert data["datasets"]["SRR1"]["exclusion"]["excluded"] is False
    assert (tmp_path / "blacklist.txt").read_text().splitlines() == ["SRR2  # test"]


def test_from_file_and_list(tmp_path, capsys):
    (tmp_path / "bad.txt").write_text("SRR5\n# comment\nSRR6\n")
    assert (
        BlacklistCommand().execute(_args(tmp_path, from_file=str(tmp_path / "bad.txt"), reason="host-dominated")) == 0
    )
    assert BlacklistCommand().execute(_args(tmp_path, list=True)) == 0
    out = capsys.readouterr().out
    assert "SRR5" in out and "host-dominated" in out


def test_registered():
    from prepare.cli.main import create_parser

    args = create_parser().parse_args(["blacklist", "--list"])
    assert args.list is True and args.blacklist_file == "blacklist.txt"
