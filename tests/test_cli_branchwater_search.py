"""Tests for the branchwater_search CLI command."""

import argparse
import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

from prepare.cli.commands.branchwater_search import BranchwaterSearchCommand
from prepare.core.constants import DEFAULT_REGISTRY_MAX_SCREENED
from prepare.core.exceptions import DataAccessError


def _args(tmp_path=None, **kwargs):
    registry_dir = Path(tmp_path) if tmp_path is not None else Path(tempfile.mkdtemp())
    base = dict(
        genome_fasta=None,
        signature=None,
        threshold=0.1,
        branchwater_folder="branchwater",
        output=None,
        server="https://s",
        registry=str(registry_dir / "prepare_registry.json"),
        registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        no_cache=False,
        refresh=False,
        max_cache_age_days=None,
    )
    base.update(kwargs)
    return argparse.Namespace(**base)


class TestBranchwaterSearchCommand:
    def test_properties(self):
        cmd = BranchwaterSearchCommand()
        assert cmd.name == "branchwater_search"
        assert "Branchwater" in cmd.help

    def test_parser_requires_one_input(self):
        parser = argparse.ArgumentParser()
        BranchwaterSearchCommand().configure_parser(parser)
        args = parser.parse_args(["--genome-fasta", "g.fna"])
        assert args.threshold == 0.1 and args.branchwater_folder == "branchwater" and args.output is None
        assert args.registry_max_screened == DEFAULT_REGISTRY_MAX_SCREENED
        assert args.no_cache is False
        assert args.refresh is False
        assert args.max_cache_age_days is None
        try:
            parser.parse_args([])
        except SystemExit as e:
            assert e.code == 2
        else:
            raise AssertionError("one of --genome-fasta/--signature must be required")

    def test_parser_accepts_cache_flags(self):
        parser = argparse.ArgumentParser()
        BranchwaterSearchCommand().configure_parser(parser)
        args = parser.parse_args(["--genome-fasta", "g.fna", "--no-cache", "--refresh", "--max-cache-age-days", "7"])
        assert args.no_cache is True
        assert args.refresh is True
        assert args.max_cache_age_days == 7

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.sketch_fasta")
    def test_genome_id_of_a_gz_fasta_drops_both_suffixes(
        self, mock_sketch, mock_search, mock_write, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        mock_sketch.return_value = {"signatures": []}
        mock_search.return_value = [("SRR1", 0.9)]
        folder = tmp_path.resolve() / "branchwater"
        assert BranchwaterSearchCommand().execute(_args(tmp_path, genome_fasta="genomes/G.fna.gz")) == 0
        mock_write.assert_called_once_with([("SRR1", 0.9)], folder / "G.csv")
        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert "G" in data["datasets"]["SRR1"]["screening"]["genomes"]

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.sketch_fasta")
    def test_fasta_default_output_path(self, mock_sketch, mock_search, mock_write, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mock_sketch.return_value = {"signatures": []}
        mock_search.return_value = [("SRR1", 0.9)]
        folder = tmp_path.resolve() / "branchwater"  # anchored on the project root (the registry's folder)
        rc = BranchwaterSearchCommand().execute(_args(tmp_path, genome_fasta="genomes/GCF_000008025.1.fna"))
        assert rc == 0
        mock_search.assert_called_once_with(
            {"signatures": []},
            0.1,
            server="https://s",
            cache_dir=folder / ".branchwater-cache",
            refresh=False,
            max_cache_age_days=None,
        )
        mock_write.assert_called_once_with([("SRR1", 0.9)], folder / "GCF_000008025.1.csv")
        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        screening = data["datasets"]["SRR1"]["screening"]
        assert screening["genomes"]["GCF_000008025.1"]["source"] == "branchwater"
        assert screening["genomes"]["GCF_000008025.1"]["containment"] == 0.9

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.load_signature")
    def test_signature_and_explicit_output(self, mock_load, mock_search, mock_write, tmp_path):
        mock_load.return_value = {"signatures": []}
        mock_search.return_value = []
        rc = BranchwaterSearchCommand().execute(_args(tmp_path, signature="wmel.sig", output=str(tmp_path / "out.csv")))
        assert rc == 0
        mock_write.assert_called_once_with([], Path(tmp_path / "out.csv"))
        assert (tmp_path / "prepare_registry.json").exists()

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index", return_value=[])
    @patch("prepare.cli.commands.branchwater_search.load_signature", return_value={"signatures": []})
    def test_zero_matches_warns(self, _load, _search, _write, caplog):
        with caplog.at_level("WARNING"):
            assert BranchwaterSearchCommand().execute(_args(signature="wmel.sig")) == 0
        assert "control genome" in caplog.text

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.load_signature")
    def test_tied_top_hit_picks_alphabetically_first_accession(self, mock_load, mock_search, _write, tmp_path, caplog):
        """A tie at the top containment must resolve the same way every run, not by list order."""
        mock_load.return_value = {"signatures": []}
        mock_search.return_value = [("SRR9", 0.9), ("SRR2", 0.9), ("SRR5", 0.9)]
        with caplog.at_level("INFO"):
            rc = BranchwaterSearchCommand().execute(
                _args(tmp_path, signature="wmel.sig", output=str(tmp_path / "out.csv"))
            )
        assert rc == 0
        assert "best containment 0.9000 (SRR2)" in caplog.text

    @patch("prepare.cli.commands.branchwater_search.load_signature", side_effect=DataAccessError("bad sig"))
    def test_error_returns_1(self, _load):
        assert BranchwaterSearchCommand().execute(_args(signature="wmel.sig")) == 1

    @patch("prepare.data.branchwater_search._session")
    @patch("prepare.cli.commands.branchwater_search.load_signature", return_value={"signatures": []})
    def test_server_ignoring_threshold_is_filtered_locally(self, _load, mock_session_factory, tmp_path, caplog):
        text = "SRA accession,containment\nSRR1,0.0009\nSRR2,0.001\nSRR3,0.0005\n"
        session = MagicMock()
        response = Mock(status_code=200, text=text)
        response.iter_lines.side_effect = lambda *args, **kwargs: iter(text.splitlines())
        session.post.return_value = response
        mock_session_factory.return_value = session
        output = tmp_path / "out.csv"
        with caplog.at_level("WARNING"):
            rc = BranchwaterSearchCommand().execute(
                _args(tmp_path, signature="wmel.sig", output=str(output), no_cache=True)
            )
        assert rc == 0
        assert "control genome" in caplog.text
        assert output.read_text().splitlines() == [
            "acc,containment,cANI,biosample,bioproject,assay_type,collection_date_sam,"
            "geo_loc_name_country_calc,organism,lat_lon"
        ]

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.load_signature")
    def test_screening_entries_are_capped(self, mock_load, mock_search, _write, tmp_path, caplog):
        """A broad search must not fill the registry; only the best matches are kept."""
        mock_load.return_value = {"signatures": []}
        mock_search.return_value = [("SRR1", 0.9), ("SRR2", 0.5), ("SRR3", 0.2)]
        with caplog.at_level("WARNING"):
            rc = BranchwaterSearchCommand().execute(
                _args(tmp_path, signature="wmel.sig", output=str(tmp_path / "out.csv"), registry_max_screened=2)
            )
        assert rc == 0
        datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
        assert sorted(datasets) == ["SRR1", "SRR2"]
        assert "2" in caplog.text

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.sketch_fasta")
    def test_no_cache_flag_disables_cache_dir(self, mock_sketch, mock_search, mock_write, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mock_sketch.return_value = {"signatures": []}
        mock_search.return_value = []
        rc = BranchwaterSearchCommand().execute(_args(tmp_path, genome_fasta="g.fna", no_cache=True))
        assert rc == 0
        _, kwargs = mock_search.call_args
        assert kwargs["cache_dir"] is None

    @patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
    @patch("prepare.cli.commands.branchwater_search.search_index")
    @patch("prepare.cli.commands.branchwater_search.sketch_fasta")
    def test_refresh_and_max_cache_age_passed_through(
        self, mock_sketch, mock_search, mock_write, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        mock_sketch.return_value = {"signatures": []}
        mock_search.return_value = []
        rc = BranchwaterSearchCommand().execute(
            _args(tmp_path, genome_fasta="g.fna", refresh=True, max_cache_age_days=3)
        )
        assert rc == 0
        _, kwargs = mock_search.call_args
        assert kwargs["refresh"] is True
        assert kwargs["max_cache_age_days"] == 3

    def test_registered(self):
        from prepare.cli.main import create_parser, register_all_commands

        register_all_commands()
        parser = create_parser()
        action = next(a for a in parser._subparsers._group_actions if getattr(a, "choices", None))
        assert "branchwater_search" in action.choices


@patch("prepare.cli.commands.branchwater_search.write_branchwater_csv")
@patch("prepare.cli.commands.branchwater_search.search_index")
@patch("prepare.cli.commands.branchwater_search.load_signature", return_value={"signatures": []})
def test_branchwater_search_drops_vanished_matches(_load, mock_search, _write, tmp_path):
    from prepare.data.registry import load_registry, query

    args = _args(tmp_path, signature=str(tmp_path / "G.sig"), output=str(tmp_path / "branchwater" / "G.csv"))
    mock_search.return_value = [("SRR1", 0.9), ("SRR2", 0.8)]
    assert BranchwaterSearchCommand().execute(args) == 0
    genome = next(iter(load_registry(args.registry).genomes))
    mock_search.return_value = [("SRR1", 0.9)]
    assert BranchwaterSearchCommand().execute(args) == 0

    assert query(load_registry(args.registry), "screened", genome) == ["SRR1"]


def _search_twice(tmp_path, first, second, threshold=0.1, seed=None):
    from prepare.data.registry import load_registry, query, registry_transaction

    args = _args(tmp_path, signature=str(tmp_path / "G.sig"), output=str(tmp_path / "branchwater" / "G.csv"))
    with (
        patch("prepare.cli.commands.branchwater_search.write_branchwater_csv"),
        patch("prepare.cli.commands.branchwater_search.load_signature", return_value={"signatures": []}),
        patch("prepare.cli.commands.branchwater_search.search_index") as mock_search,
    ):
        args.threshold = 0.0
        mock_search.return_value = first
        assert BranchwaterSearchCommand().execute(args) == 0
        genome = next(iter(load_registry(args.registry).genomes))
        if seed is not None:
            with registry_transaction(args.registry) as registry:
                seed(registry, genome)
        args.threshold = threshold
        mock_search.return_value = second
        assert BranchwaterSearchCommand().execute(args) == 0
    return sorted(query(load_registry(args.registry), "screened", genome))


def test_an_empty_branchwater_reply_keeps_the_entries(tmp_path):
    """Review: the public index has answered nothing for known positives; that must not erase a genome."""
    assert _search_twice(tmp_path, [("SRR1", 0.9), ("SRR2", 0.5)], []) == ["SRR1", "SRR2"]


def test_a_higher_threshold_search_keeps_entries_below_it(tmp_path):
    """Review: a search at 0.1 says nothing about accessions below 0.1."""
    assert _search_twice(tmp_path, [("SRR1", 0.9), ("SRR2", 0.05)], [("SRR1", 0.9)]) == ["SRR1", "SRR2"]


def test_a_branchwater_search_leaves_entries_from_match_tables_alone(tmp_path):
    from prepare.data.registry import record_screening

    def seed(registry, genome):
        record_screening(registry, "SRR3", genome, 0.4, None, "matches", 0.0, None)

    assert _search_twice(tmp_path, [("SRR1", 0.9)], [("SRR1", 0.9)], seed=seed) == ["SRR1", "SRR3"]


def _genomes_after_search(tmp_path, output):
    from prepare.data.registry import load_registry

    args = _args(tmp_path, signature=str(tmp_path / "GCF_1.sig"), output=output)
    with (
        patch("prepare.cli.commands.branchwater_search.write_branchwater_csv"),
        patch("prepare.cli.commands.branchwater_search.load_signature", return_value={"signatures": []}),
        patch("prepare.cli.commands.branchwater_search.search_index", return_value=[("SRR1", 0.9)]),
    ):
        assert BranchwaterSearchCommand().execute(args) == 0
    return list(load_registry(args.registry).genomes)


def test_a_custom_output_records_under_its_stem(tmp_path):
    """Review: parse_containment names the column after the CSV stem, so the registry must too."""
    assert _genomes_after_search(tmp_path, str(tmp_path / "branchwater" / "foo.csv")) == ["foo"]


def test_default_output_records_under_the_genome_name(tmp_path):
    assert _genomes_after_search(tmp_path, None) == ["GCF_1"]
