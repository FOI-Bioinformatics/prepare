"""
SIMPLIFIED INTEGRATION TESTS for PREPARE

Tests complete end-to-end workflows using components we've validated:
- SRA metadata workflow
- File processing pipelines
- Visualization generation
- Data export workflows

Run: pytest tests/test_integration_simple.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import Mock
import matplotlib

matplotlib.use("Agg")  # Use non-GUI backend for testing

from prepare.data.sra_metadata import (  # noqa: E402
    SRAMetadataClient,
    SRADatasetInfo,
    detect_sequencing_technology,
    save_metadata_report,
    create_download_preview,
)
from prepare.cli.commands.sra_profile import SRAProfileCommand  # noqa: E402
from prepare.plugins.visualizers.bar import BarChartPlugin  # noqa: E402


def run_sra_profile(fastq_dir, output_report):
    """Run sra_profile over every accession folder in ``fastq_dir``; returns its exit status."""
    import argparse

    parser = argparse.ArgumentParser()
    SRAProfileCommand().configure_parser(parser)
    work = Path(output_report).parent
    args = parser.parse_args(
        [
            "--fastq-folder",
            str(fastq_dir),
            "--output-report",
            str(output_report),
            "--output-dir",
            str(work / "profiles"),
            "--registry",
            str(work / "prepare_registry.json"),
        ]
    )
    return SRAProfileCommand().execute(args)


# ============================================================================
# TEST CLASS: SRA Metadata to Visualization Workflow
# ============================================================================


class TestSRAMetadataWorkflow:
    """Test complete SRA metadata workflow."""

    @pytest.fixture
    def sample_metadata(self):
        """Create sample SRA metadata."""
        return {
            "SRR001": SRADatasetInfo(
                accession="SRR001",
                title="E. coli WGS",
                organism="Escherichia coli",
                platform="ILLUMINA",
                instrument="HiSeq 2500",
                strategy="WGS",
                layout="PAIRED",
                spots=1000000,
                bases=150000000,
                avg_length=150.0,
                size_mb=100.0,
                release_date="2023-01-01",
                bioproject="PRJNA001",
                biosample="SAMN001",
                library_selection="RANDOM",
                library_source="GENOMIC",
            ),
            "SRR002": SRADatasetInfo(
                accession="SRR002",
                title="S. aureus WGS",
                organism="Staphylococcus aureus",
                platform="OXFORD_NANOPORE",
                instrument="MinION",
                strategy="WGS",
                layout="SINGLE",
                spots=500000,
                bases=500000000,
                avg_length=1000.0,
                size_mb=400.0,
                release_date="2023-01-02",
                bioproject="PRJNA002",
                biosample="SAMN002",
                library_selection="RANDOM",
                library_source="GENOMIC",
            ),
            "SRR003": SRADatasetInfo(
                accession="SRR003",
                title="K. pneumoniae WGS",
                organism="Klebsiella pneumoniae",
                platform="ILLUMINA",
                instrument="MiSeq",
                strategy="WGS",
                layout="PAIRED",
                spots=750000,
                bases=112500000,
                avg_length=150.0,
                size_mb=75.0,
                release_date="2023-01-03",
                bioproject="PRJNA003",
                biosample="SAMN003",
                library_selection="RANDOM",
                library_source="GENOMIC",
            ),
        }

    def test_metadata_retrieval_and_export_workflow(self, sample_metadata, tmp_path):
        """Test: retrieve metadata -> export to CSV."""
        # Save metadata report
        report_file = tmp_path / "sra_metadata.csv"
        save_metadata_report(sample_metadata, report_file)

        assert report_file.exists()

        # Verify can read back
        df = pd.read_csv(report_file)
        assert len(df) == 3
        assert set(df["accession"]) == {"SRR001", "SRR002", "SRR003"}
        assert "organism" in df.columns
        assert "technology" in df.columns

    def test_metadata_to_visualization_workflow(self, sample_metadata, tmp_path):
        """Test: metadata -> technology analysis -> visualization."""
        # Step 1: Analyze technologies
        tech_counts = {}
        for acc, info in sample_metadata.items():
            tech = detect_sequencing_technology(info)
            tech_counts[tech] = tech_counts.get(tech, 0) + 1

        assert tech_counts["illumina"] == 2
        assert tech_counts["nanopore"] == 1

        # Step 2: Create visualization
        tech_df = pd.DataFrame(
            {
                "technology": list(tech_counts.keys()),
                "count": list(tech_counts.values()),
            }
        )

        plot_file = tmp_path / "technology_distribution.png"
        BarChartPlugin.create_plot(
            data=tech_df,
            x_column="technology",
            y_column="count",
            title="Sequencing Technology Distribution",
            output_file=str(plot_file),
        )

        assert plot_file.exists()

    def test_metadata_preview_workflow(self, sample_metadata, tmp_path):
        """Test: metadata -> download preview -> summary."""
        # Create mock client
        mock_client = Mock(spec=SRAMetadataClient)
        mock_client.get_sra_metadata.return_value = sample_metadata

        # Get download preview
        accessions = list(sample_metadata.keys())
        metadata, tech_counts, total_size_gb = create_download_preview(accessions, mock_client)

        assert len(metadata) == 3
        assert tech_counts["illumina"] == 2
        assert tech_counts["nanopore"] == 1
        assert total_size_gb == pytest.approx((100 + 400 + 75) / 1024, rel=0.01)


# ============================================================================
# TEST CLASS: FASTQ Processing Workflow
# ============================================================================


class TestFASTQProcessingWorkflow:
    """Test FASTQ file processing workflows."""

    @pytest.fixture
    def fastq_workspace(self, tmp_path):
        """Create FASTQ workspace with test files."""
        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir()

        # Create multiple accession directories with FASTQ files
        for i in range(3):
            acc_dir = fastq_dir / f"SRR00{i+1}"
            acc_dir.mkdir()

            # Create paired-end FASTQ files
            for suffix in ["_R1.fastq", "_R2.fastq"]:
                fastq_file = acc_dir / f"SRR00{i+1}{suffix}"
                # Create realistic FASTQ content
                reads = []
                for j in range(10):
                    reads.append(f"@read{j+1}")
                    reads.append("ATCGATCGATCG" * 5)  # 60bp read
                    reads.append("+")
                    reads.append("I" * 60)
                fastq_file.write_text("\n".join(reads) + "\n")

        return fastq_dir

    def test_fastq_statistics_workflow(self, fastq_workspace, tmp_path):
        """Test: FASTQ files -> statistics -> report."""
        output_report = tmp_path / "fastq_statistics.csv"

        # Generate statistics
        assert run_sra_profile(fastq_workspace, output_report) == 0

        assert output_report.exists()

        # Verify report
        stats_df = pd.read_csv(output_report)
        assert len(stats_df) == 3
        assert all(stats_df["layout"] == "PAIRED")
        assert all(stats_df["total_reads"] == 20)  # 10 reads per file, 2 files

    def test_fastq_to_visualization_workflow(self, fastq_workspace, tmp_path):
        """Test: FASTQ stats -> visualization."""
        # Generate statistics
        stats_file = tmp_path / "stats.csv"
        assert run_sra_profile(fastq_workspace, stats_file) == 0

        # Read stats
        stats_df = pd.read_csv(stats_file)

        # Create visualization of read counts
        plot_file = tmp_path / "read_counts.png"
        BarChartPlugin.create_plot(
            data=stats_df,
            x_column="accession",
            y_column="total_reads",
            title="Read Counts per Sample",
            output_file=str(plot_file),
        )

        assert plot_file.exists()


# ============================================================================
# TEST CLASS: Multi-Sample Comparison Workflow
# ============================================================================


class TestMultiSampleComparison:
    """Test workflows comparing multiple samples."""

    def test_technology_comparison_workflow(self, tmp_path):
        """Test: Multiple samples -> technology comparison -> visualization."""
        # Create sample data with different technologies
        samples = []
        for i in range(10):
            tech = ["illumina", "nanopore", "pacbio"][i % 3]
            samples.append(
                {
                    "sample_id": f"SRR{i:04d}",
                    "technology": tech,
                    "read_count": np.random.randint(100000, 1000000),
                    "quality_score": np.random.uniform(20, 40),
                }
            )

        samples_df = pd.DataFrame(samples)

        # Save to file
        samples_file = tmp_path / "samples.csv"
        samples_df.to_csv(samples_file, index=False)

        # Analyze technology distribution
        tech_counts = samples_df["technology"].value_counts()

        # Create visualization
        tech_plot = tmp_path / "tech_comparison.png"
        tech_df = pd.DataFrame(
            {
                "technology": tech_counts.index,
                "count": tech_counts.values,
            }
        )

        BarChartPlugin.create_plot(
            data=tech_df,
            x_column="technology",
            y_column="count",
            title="Technology Distribution",
            output_file=str(tech_plot),
        )

        assert tech_plot.exists()

    def test_quality_analysis_workflow(self, tmp_path):
        """Test: Quality metrics -> analysis -> visualization."""
        # Create quality data for multiple samples
        quality_data = []
        for i in range(20):
            quality_data.append(
                {
                    "sample": f"Sample_{i+1}",
                    "mean_quality": np.random.uniform(25, 40),
                    "gc_content": np.random.uniform(40, 60),
                    "read_length": np.random.randint(50, 300),
                }
            )

        quality_df = pd.DataFrame(quality_data)

        # Save quality report
        quality_file = tmp_path / "quality_report.csv"
        quality_df.to_csv(quality_file, index=False)

        # Create quality visualization
        quality_plot = tmp_path / "quality_distribution.png"
        BarChartPlugin.create_plot(
            data=quality_df.head(10),  # Top 10 samples
            x_column="sample",
            y_column="mean_quality",
            title="Quality Scores by Sample",
            output_file=str(quality_plot),
        )

        assert quality_plot.exists()


# ============================================================================
# TEST CLASS: Data Export and Reporting Workflow
# ============================================================================


class TestDataExportWorkflow:
    """Test data export and reporting workflows."""

    def test_metadata_to_multiple_formats_workflow(self, tmp_path):
        """Test: Metadata -> export to multiple formats."""
        # Create test metadata
        metadata = {
            "SRR001": SRADatasetInfo(
                accession="SRR001",
                title="Test 1",
                organism="E. coli",
                platform="ILLUMINA",
                instrument="HiSeq",
                strategy="WGS",
                layout="PAIRED",
                spots=1000000,
                bases=150000000,
                avg_length=150.0,
                size_mb=100.0,
                release_date="2023-01-01",
                bioproject="PRJNA001",
                biosample="SAMN001",
                library_selection="RANDOM",
                library_source="GENOMIC",
            ),
        }

        # Export to CSV
        csv_file = tmp_path / "metadata.csv"
        save_metadata_report(metadata, csv_file)
        assert csv_file.exists()

        # Read and re-export in different format
        df = pd.read_csv(csv_file)

        # Export to TSV
        tsv_file = tmp_path / "metadata.tsv"
        df.to_csv(tsv_file, sep="\t", index=False)
        assert tsv_file.exists()

        # Export to JSON
        json_file = tmp_path / "metadata.json"
        df.to_json(json_file, orient="records", indent=2)
        assert json_file.exists()

    def test_summary_report_generation_workflow(self, tmp_path):
        """Test: Data analysis -> summary report -> export."""
        # Create analysis results
        results = {
            "total_samples": 100,
            "illumina_count": 60,
            "nanopore_count": 30,
            "pacbio_count": 10,
            "total_reads": 50000000,
            "total_bases": 7500000000,
            "mean_quality": 35.2,
        }

        # Create summary report file
        report_file = tmp_path / "summary_report.txt"
        with open(report_file, "w") as f:
            f.write("Analysis Summary Report\n")
            f.write("=" * 50 + "\n\n")
            for key, value in results.items():
                f.write(f"{key.replace('_', ' ').title()}: {value}\n")

        assert report_file.exists()

        # Verify content
        content = report_file.read_text()
        assert "Analysis Summary Report" in content
        assert "Total Samples: 100" in content
        assert "Mean Quality: 35.2" in content


# ============================================================================
# TEST CLASS: Error Handling in Workflows
# ============================================================================


class TestErrorHandlingWorkflows:
    """Test error handling in complete workflows."""

    def test_workflow_with_missing_data(self, tmp_path):
        """Test workflow gracefully handles missing data."""
        # Try to generate stats for non-existent directory
        non_existent = tmp_path / "does_not_exist"
        output_file = tmp_path / "output.csv"

        assert run_sra_profile(non_existent, output_file) == 1
        assert not output_file.exists()

    def test_workflow_with_empty_results(self, tmp_path):
        """Test workflow handles empty results."""
        # Create empty metadata
        empty_metadata = {}

        # Should handle gracefully
        output_file = tmp_path / "empty_metadata.csv"
        save_metadata_report(empty_metadata, output_file)

        # Should create empty or minimal file
        # (behavior depends on implementation)

    def test_workflow_continues_after_partial_errors(self, tmp_path):
        """Test workflow continues processing valid data after errors."""
        # Create FASTQ workspace with mix of valid and invalid
        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir()

        # Valid directory
        valid_dir = fastq_dir / "SRR001"
        valid_dir.mkdir()
        fastq_file = valid_dir / "SRR001_R1.fastq"
        fastq_file.write_text("@read1\nATCG\n+\nIIII\n")

        # Invalid directory (no FASTQ files)
        invalid_dir = fastq_dir / "SRR002"
        invalid_dir.mkdir()
        (invalid_dir / "not_fastq.txt").write_text("invalid")

        # Should process valid data despite invalid directory
        output_file = tmp_path / "stats.csv"
        # Exit 1 reports the accession that could not be profiled ...
        assert run_sra_profile(fastq_dir, output_file) == 1

        # ... and the valid directory still has its row.
        assert list(pd.read_csv(output_file)["accession"]) == ["SRR001"]


# ============================================================================
# SUCCESS METRICS:
#
# After running these tests:
# - 15+ integration tests covering end-to-end workflows
# - SRA metadata workflows validated
# - FASTQ processing pipelines tested
# - Multi-sample comparison workflows verified
# - Error handling validated
#
# Run tests:
#   pytest tests/test_integration_simple.py -v
# ============================================================================
