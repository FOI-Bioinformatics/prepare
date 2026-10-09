"""
Unit tests for PREPARE data models.
"""

from prepare.core.models import (
    Containment,
    GenomeInfo,
    SRAMetadata,
    ContainmentSummary,
)


def test_containment_model():
    """Test the Containment data model."""
    # Create a basic containment object
    containment = Containment(accession="SRR123456", value=0.95, genome_id="GCF_123456.1")

    assert containment.accession == "SRR123456"
    assert containment.value == 0.95
    assert containment.genome_id == "GCF_123456.1"
    assert containment.additional_data == {}

    # Create containment with additional data
    containment_with_data = Containment(
        accession="SRR789012",
        value=0.87,
        genome_id="GCF_789012.1",
        additional_data={"biosample": "SAMN123456", "organism": "Escherichia coli"},
    )

    assert containment_with_data.accession == "SRR789012"
    assert containment_with_data.value == 0.87
    assert containment_with_data.genome_id == "GCF_789012.1"
    assert containment_with_data.additional_data["biosample"] == "SAMN123456"
    assert containment_with_data.additional_data["organism"] == "Escherichia coli"


def test_genome_info_model():
    """Test the GenomeInfo data model."""
    # Create a basic genome info object
    genome = GenomeInfo(genome_id="GCF_123456.1")

    assert genome.genome_id == "GCF_123456.1"
    assert genome.path is None
    assert genome.description is None
    assert genome.metadata == {}

    # Create genome info with all fields
    from pathlib import Path
    import datetime

    genome_full = GenomeInfo(
        genome_id="GCF_789012.1",
        path=Path("/genomes/GCF_789012.1.fasta"),
        description="Test genome",
        source="NCBI",
        date_added=datetime.datetime.now(),
        metadata={"organism": "Escherichia coli", "assembly_level": "Complete Genome"},
    )

    assert genome_full.genome_id == "GCF_789012.1"
    assert genome_full.path.name == "GCF_789012.1.fasta"
    assert genome_full.description == "Test genome"
    assert genome_full.source == "NCBI"
    assert genome_full.date_added is not None
    assert genome_full.metadata["organism"] == "Escherichia coli"


def test_sra_metadata_model():
    """Test the SRAMetadata data model."""
    # Create a basic SRA metadata object
    metadata = SRAMetadata(accession="SRR123456")

    assert metadata.accession == "SRR123456"
    assert metadata.biosample is None
    assert metadata.bioproject is None
    assert metadata.organism is None
    assert metadata.attributes == {}

    # Create SRA metadata with all fields
    metadata_full = SRAMetadata(
        accession="SRR789012",
        biosample="SAMN789012",
        bioproject="PRJNA789012",
        organism="Escherichia coli",
        collection_date="2020-01-01",
        location="USA: California",
        latitude_longitude="37.7749 N 122.4194 W",
        assay_type="WGS",
        attributes={"host": "Homo sapiens", "isolation_source": "stool"},
    )

    assert metadata_full.accession == "SRR789012"
    assert metadata_full.biosample == "SAMN789012"
    assert metadata_full.bioproject == "PRJNA789012"
    assert metadata_full.organism == "Escherichia coli"
    assert metadata_full.collection_date == "2020-01-01"
    assert metadata_full.location == "USA: California"
    assert metadata_full.latitude_longitude == "37.7749 N 122.4194 W"
    assert metadata_full.assay_type == "WGS"
    assert metadata_full.attributes["host"] == "Homo sapiens"

    # Test the country property
    assert metadata_full.country == "USA"

    # Test with different location format
    metadata_different_format = SRAMetadata(accession="SRR456789", location="France")
    assert metadata_different_format.country == "France"


def test_containment_summary_model():
    """Test the ContainmentSummary data model."""
    # Create a basic containment summary object
    summary = ContainmentSummary()

    assert summary.thresholds == []
    assert summary.counts == []
    assert summary.max_containment == {}
    assert summary.genome_to_samples == {}
    assert summary.sample_to_genomes == {}

    # Create containment summary with data
    summary_with_data = ContainmentSummary(
        thresholds=[0.9, 0.8, 0.7],
        counts=[10, 20, 30],
        max_containment={"SRR123456": 0.95, "SRR789012": 0.85},
        genome_to_samples={
            "GCF_123456.1": ["SRR123456", "SRR789012"],
            "GCF_789012.1": ["SRR123456"],
        },
        sample_to_genomes={
            "SRR123456": ["GCF_123456.1", "GCF_789012.1"],
            "SRR789012": ["GCF_123456.1"],
        },
    )

    assert summary_with_data.thresholds == [0.9, 0.8, 0.7]
    assert summary_with_data.counts == [10, 20, 30]
    assert summary_with_data.max_containment["SRR123456"] == 0.95
    assert "GCF_123456.1" in summary_with_data.genome_to_samples
    assert "SRR123456" in summary_with_data.sample_to_genomes
    assert len(summary_with_data.genome_to_samples["GCF_123456.1"]) == 2
    assert len(summary_with_data.sample_to_genomes["SRR123456"]) == 2
