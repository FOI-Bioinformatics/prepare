"""Optional dependencies: the CLI loads without them and a missing one names its extra."""

import argparse
import pathlib
import sys
import tomllib

import pytest

from prepare.core.exceptions import ConfigurationError
from prepare.core.optional import require

OPTIONAL_MODULES = (
    "sklearn",
    "scipy",
    "plotly",
    "jinja2",
    "seaborn",
    "cartopy",
    "sourmash",
    "umap",
    "statsmodels",
    "networkx",
    "upsetplot",
)

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def _args(tmp_path):
    table = tmp_path / "abundance.csv"
    table.write_text("sample,sp1,sp2\ns1,10,0\ns2,3,7\n")
    return argparse.Namespace(
        abundance_file=str(table),
        metadata_file=None,
        output_dir=str(tmp_path / "out"),
        alpha_metrics=["shannon"],
        beta_metric="bray_curtis",
        permanova_formula=None,
    )


@pytest.fixture
def fresh_prepare_modules(monkeypatch):
    """Drop every prepare module so the next import runs module-level code again.

    monkeypatch restores the modules it removed; the copies imported during the test are
    removed here, so later tests keep the original classes (for pytest.raises identity).
    """
    for mod in list(sys.modules):
        if mod.startswith("prepare"):
            monkeypatch.delitem(sys.modules, mod)
    yield
    for mod in list(sys.modules):
        if mod.startswith("prepare"):
            del sys.modules[mod]


def test_cli_imports_without_optional_packages(monkeypatch, fresh_prepare_modules):
    for name in OPTIONAL_MODULES:
        monkeypatch.setitem(sys.modules, name, None)
    import prepare.cli.main as main  # must not raise

    parser = main.create_parser()
    assert parser is not None
    args = parser.parse_args(["diversity_analysis", "--abundance-file", "x.csv"])
    assert args.abundance_file == "x.csv"
    assert parser.parse_args(["sra_sketch", "--accession", "SRR1"]).accession == ["SRR1"]


def test_diversity_analysis_names_the_extra_when_sklearn_is_missing(monkeypatch, tmp_path):
    from prepare.cli.commands.advanced_analysis import DiversityAnalysisCommand

    monkeypatch.setitem(sys.modules, "sklearn", None)
    with pytest.raises(ConfigurationError) as exc:
        DiversityAnalysisCommand().execute(_args(tmp_path))
    assert "prepare[analysis]" in str(exc.value) and sys.executable in str(exc.value)
    # The message names the distribution pip knows, not the import name.
    assert "needs the 'scikit-learn' package" in str(exc.value)
    # The check runs before any work, so no partial output is left behind.
    assert not (tmp_path / "out" / "alpha_diversity.csv").exists()


def test_pyproject_has_no_dead_dependencies():
    text = (REPO_ROOT / "pyproject.toml").read_text()
    for dead in ("statsmodels", "umap-learn", "networkx", "upsetplot", "seaborn"):
        assert dead not in text


def test_core_dependencies_are_the_core_set():
    project = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())["project"]
    names = {dep.split(">")[0].split("=")[0].strip() for dep in project["dependencies"]}
    assert names == {"pandas", "numpy", "matplotlib", "biopython", "lxml", "requests"}
    for extra in ("analysis", "interactive", "maps", "sourmash", "all"):
        assert extra in project["optional-dependencies"]


def test_require_returns_the_module_when_present():
    plotly = pytest.importorskip("plotly")
    assert require("plotly", "interactive", "x") is plotly


def test_require_names_extra_and_interpreter_when_missing(monkeypatch):
    monkeypatch.setitem(sys.modules, "plotly", None)
    with pytest.raises(ConfigurationError) as exc:
        require("plotly.graph_objects", "interactive", "An interactive plot")
    message = str(exc.value)
    assert message.startswith("An interactive plot needs the 'plotly' package.")
    assert f"{sys.executable} -m pip install 'prepare[interactive]'" in message


def test_sra_dashboard_without_plotly_is_an_error_not_a_degraded_file(monkeypatch, tmp_path):
    from prepare.sra.reporting import SRAReportGenerator

    monkeypatch.setitem(sys.modules, "plotly", None)
    generator = SRAReportGenerator(tmp_path / "reports")
    with pytest.raises(ConfigurationError, match=r"prepare\[interactive\]"):
        generator.generate_quality_dashboard(["SRR000001"])
    assert list((tmp_path / "reports").iterdir()) == []


def test_explorer_without_plotly_is_an_error(monkeypatch, tmp_path):
    from prepare.visualization import explorer

    monkeypatch.setitem(sys.modules, "plotly", None)
    with pytest.raises(ConfigurationError, match=r"prepare\[interactive\]"):
        explorer.require_explorer_packages()
