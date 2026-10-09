"""Tests for prepare/visualization/explorer.py.

Covers HTML generation, empty-data handling, metadata integration,
min_containment filtering, and output validation.

Run: pytest tests/test_explorer.py -v
"""

import pytest
import pandas as pd

for _module in ("plotly", "jinja2"):
    pytest.importorskip(_module)

from prepare.core.models import TaxonomyInfo  # noqa: E402
from prepare.visualization.explorer import (  # noqa: E402
    generate_containment_explorer,
    _build_long_dataframe,
    _build_summary_data,
    _build_table_html,
)


@pytest.fixture()
def containment_df():
    return pd.DataFrame(
        {"genA": [0.5, 0.1, 0.8], "genB": [0.2, 0.9, 0.0]},
        index=["SRR001", "SRR002", "SRR003"],
    )


@pytest.fixture()
def taxonomy():
    return {
        "genA": TaxonomyInfo(genome_id="genA", family="Fam1", genus="Gen1", species="Sp1"),
        "genB": TaxonomyInfo(genome_id="genB", family="Fam2", genus="Gen2", species="Sp2"),
    }


@pytest.fixture()
def metadata_df():
    return pd.DataFrame(
        {"location": ["USA", "Sweden", "Japan"], "collection_date": ["2024-01", "2024-02", "2024-03"]},
        index=["SRR001", "SRR002", "SRR003"],
    )


class TestGenerateContainmentExplorer:

    def test_generates_html_file(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "explorer.html"
        result = generate_containment_explorer(containment_df, taxonomy, output_file=out)
        assert result.exists()
        assert result.suffix == ".html"

    def test_html_contains_expected_sections(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "test.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out)
        html = out.read_text()
        assert "resultsTable" in html
        assert "Taxonomy overview" in html
        assert "plotly-graph-div" in html  # charts rendered
        assert "mq-heatbar" in html  # containment heat-bar cells
        assert "SRR001" in html
        assert "filterTable" in html

    def test_with_metadata(self, containment_df, taxonomy, metadata_df, tmp_path):
        out = tmp_path / "meta.html"
        generate_containment_explorer(containment_df, taxonomy, metadata=metadata_df, output_file=out)
        html = out.read_text()
        assert "USA" in html
        assert "Sweden" in html

    def test_min_containment_filtering(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "filtered.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out, min_containment=0.5)
        html = out.read_text()
        # SRR002 has genA=0.1 and genB=0.9. genA should be excluded but genB included.
        assert "SRR002" in html  # still present via genB=0.9

    def test_empty_taxonomy(self, containment_df, tmp_path):
        out = tmp_path / "empty_tax.html"
        generate_containment_explorer(containment_df, {}, output_file=out)
        html = out.read_text()
        assert "Unknown" in html
        assert "resultsTable" in html

    def test_custom_title(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "titled.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out, title="Custom Title")
        html = out.read_text()
        assert "Custom Title" in html

    def test_output_in_new_subdirectory(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "sub" / "dir" / "out.html"
        result = generate_containment_explorer(containment_df, taxonomy, output_file=out)
        assert result.exists()


class TestBuildLongDataframe:

    def test_basic_conversion(self, containment_df, taxonomy):
        df = _build_long_dataframe(containment_df, taxonomy, None, 0.0)
        assert "sample" in df.columns
        assert "genome" in df.columns
        assert "family" in df.columns
        # Only rows with containment > 0
        assert len(df) == 5  # 6 total cells, one is 0.0

    def test_min_containment(self, containment_df, taxonomy):
        df = _build_long_dataframe(containment_df, taxonomy, None, 0.5)
        assert all(df["containment"] >= 0.5)

    def test_min_containment_is_inclusive(self, containment_df, taxonomy):
        # genA at SRR001 is exactly 0.5 and must be kept when min_containment=0.5.
        df = _build_long_dataframe(containment_df, taxonomy, None, 0.5)
        assert ((df["sample"] == "SRR001") & (df["genome"] == "genA")).any()


class TestBuildSummaryData:

    def test_summary_values(self, containment_df, taxonomy):
        long_df = _build_long_dataframe(containment_df, taxonomy, None, 0.0)
        summary = _build_summary_data(long_df, containment_df, taxonomy)
        assert summary["total_samples"] == 3
        assert summary["total_genomes"] == 2
        assert summary["num_families"] == 2
        assert summary["containment_max"] == 0.9

    def test_empty_long_df(self, containment_df, taxonomy):
        empty = pd.DataFrame()
        summary = _build_summary_data(empty, containment_df, taxonomy)
        assert summary["total_samples"] == 0


class TestBuildTableHtml:

    def test_empty_df(self):
        assert _build_table_html(pd.DataFrame()) == ""

    def test_rows_generated(self, containment_df, taxonomy):
        long_df = _build_long_dataframe(containment_df, taxonomy, None, 0.0)
        html = _build_table_html(long_df)
        assert "<tr>" in html
        assert "SRR001" in html
