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
    _sunburst_figure,
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
        assert "pre-heatbar" in html  # containment heat-bar cells
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

    def test_no_taxonomy_hides_the_taxonomy_sections(self, containment_df, tmp_path):
        """G6: without a taxonomy map the page says so instead of charting 'Unknown'."""
        out = tmp_path / "empty_tax.html"
        generate_containment_explorer(containment_df, {}, output_file=out)
        html = out.read_text()
        assert "resultsTable" in html
        assert "Taxonomy overview" not in html and "Distributions" not in html
        assert '<p class="k">Families</p>' not in html and '<p class="k">Genera</p>' not in html
        assert "Top families" not in html
        assert "enrich_taxonomy" in html

    def test_an_empty_result_explains_itself(self, containment_df, taxonomy, tmp_path):
        """G10: no row above --min-containment: the page says so and still counts the genome columns."""
        out = tmp_path / "empty.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out, min_containment=0.95)
        html = out.read_text()
        import re

        assert re.search(r"No sample reaches containment 0\.950;\s+lower --min-containment", html)
        assert "<b>2</b> genomes" in html

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
        assert summary["total_genomes"] == 2  # the genome columns are still there (G10)


class TestBuildTableHtml:

    def test_empty_df(self):
        assert _build_table_html(pd.DataFrame()) == ""

    def test_rows_generated(self, containment_df, taxonomy):
        long_df = _build_long_dataframe(containment_df, taxonomy, None, 0.0)
        html = _build_table_html(long_df)
        assert "<tr>" in html
        assert "SRR001" in html


class TestExplorerEscapingAndMetadata:

    def test_family_name_is_escaped(self, containment_df, tmp_path):
        evil = "<img src=x onerror=alert(2)>"
        taxonomy = {
            "genA": TaxonomyInfo(genome_id="genA", family=evil, genus="Gen1", species="Sp1"),
            "genB": TaxonomyInfo(genome_id="genB", family="Fam2", genus="Gen2", species="Sp2"),
        }
        out = tmp_path / "x.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out)
        html = out.read_text()
        assert "&lt;img" in html
        assert "<img" not in html

    def test_title_is_escaped_once(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "t.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out, title="A & <b>")
        html = out.read_text()
        assert "A &amp; &lt;b&gt;" in html
        assert "&amp;amp;" not in html

    def test_prepare_metadata_columns_fill_location_and_date(self, containment_df, taxonomy):
        meta = pd.DataFrame(
            {"geo_loc_name_country_calc": ["Peru", "Kenya", "Fiji"], "collection_date_sam": ["2020", "2021", "2022"]},
            index=["SRR001", "SRR002", "SRR003"],
        )
        long_df = _build_long_dataframe(containment_df, taxonomy, meta, 0.0)
        row = long_df[long_df["sample"] == "SRR001"].iloc[0]
        assert row["location"] == "Peru"
        assert row["date"] == "2020"

    def test_nan_metadata_is_an_empty_cell(self, containment_df, taxonomy):
        meta = pd.DataFrame(
            {"location": [float("nan"), "Sweden", "Japan"], "collection_date": ["2024-01", None, "2024-03"]},
            index=["SRR001", "SRR002", "SRR003"],
        )
        long_df = _build_long_dataframe(containment_df, taxonomy, meta, 0.0)
        assert "nan" not in set(long_df["location"]) | set(long_df["date"])
        assert long_df[long_df["sample"] == "SRR001"].iloc[0]["location"] == ""

    def test_duplicated_metadata_index_does_not_raise(self, containment_df, taxonomy):
        meta = pd.DataFrame(
            {"location": ["USA", "Peru", "Japan"], "collection_date": ["2024-01", "2024-02", "2024-03"]},
            index=["SRR001", "SRR001", "SRR003"],
        )
        long_df = _build_long_dataframe(containment_df, taxonomy, meta, 0.0)
        assert long_df[long_df["sample"] == "SRR001"].iloc[0]["location"] == "USA"


class TestExplorerScaling:

    def test_large_table_builds_quickly_and_is_capped(self, tmp_path):
        import time

        import numpy as np

        rng = np.random.default_rng(0)
        genomes = [f"g{i}" for i in range(50)]
        wide = pd.DataFrame(rng.random((300, 50)) + 0.01, index=[f"S{i}" for i in range(300)], columns=genomes)
        tax = {
            g: TaxonomyInfo(genome_id=g, family=f"F{i % 5}", genus=f"G{i % 10}", species=f"s{i}")
            for i, g in enumerate(genomes)
        }
        start = time.perf_counter()
        long_df = _build_long_dataframe(wide, tax, None, 0.0)
        assert time.perf_counter() - start < 2.0
        assert len(long_df) == 300 * 50
        assert _build_table_html(long_df).count("<tr>") == 5000

    def test_header_note_only_when_truncated(self, containment_df, taxonomy, tmp_path):
        out = tmp_path / "n.html"
        generate_containment_explorer(containment_df, taxonomy, output_file=out)
        assert "showing top" not in out.read_text()

    def test_truncated_table_has_note(self, tmp_path, monkeypatch):
        import prepare.visualization.explorer as ex

        monkeypatch.setattr(ex, "EXPLORER_TABLE_ROWS", 3)
        wide = pd.DataFrame({"gA": [0.5, 0.4, 0.3, 0.2]}, index=["a", "b", "c", "d"])
        tax = {"gA": TaxonomyInfo(genome_id="gA", family="F", genus="G", species="S")}
        out = tmp_path / "t.html"
        ex.generate_containment_explorer(wide, tax, output_file=out)
        html = out.read_text()
        assert "showing top 3 of 4 rows" in html
        assert html.count("<tr>") == 3 + 1  # header row plus the capped body rows

    def test_heatmap_height_is_capped(self):
        from prepare.visualization.explorer import _build_heatmap

        samples = [f"S{i}" for i in range(400)]
        long_df = pd.DataFrame(
            {"sample": samples, "family": "F", "containment": 0.5, "genus": "G", "species": "s", "genome": "g"}
        )
        assert '"height":2000' in _build_heatmap(long_df).replace(" ", "")


def test_sunburst_counts_hits_so_the_hierarchy_adds_up():
    """A sample holding two genomes of one family is two hits, and the chart says hits, not samples."""
    long_df = pd.DataFrame(
        {
            "sample": ["S1", "S1", "S2", "S2"],
            "genome": ["g1", "g2", "g1", "g2"],
            "containment": [0.5, 0.7, 0.6, 0.8],
            "family": ["F"] * 4,
            "genus": ["G"] * 4,
            "species": ["s1", "s2", "s1", "s2"],
        }
    )
    fig = _sunburst_figure(long_df)
    trace = fig.data[0]
    values = dict(zip(trace.labels, trace.values))
    assert values["F"] == 4 and values["s1"] == 2
    assert "hits" in fig.layout.title.text.lower()
    assert "sample_count" not in fig.to_json()
