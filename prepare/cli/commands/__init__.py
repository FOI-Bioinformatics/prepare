"""
CLI commands package.

This package contains modular CLI command implementations.
"""

from .branchwater import UseBranchwaterCommand, ExtractBranchwaterMetadataCommand
from .containment import ParseContainmentCommand, PlotContainmentCommand
from .metadata import (
    DownloadMetadataCommand,
    ParseMetadataCommand,
    CheckMetadataAttributesCommand,
    CountMetadataCommand,
    PlotMetadataCountsCommand,
)
from .samples import SingleSampleCommand
from .sra import DownloadSraCommand
from .status import StatusCommand
from .results import ResultsTableCommand
from .read_extraction import ExtractTargetReadsCommand
from .test_data import DownloadTestGenomeCommand
from .genome import GenomeSearchCommand, GenomeDownloadCommand, GenomePrepareCommand
from .explore import (
    EnrichTaxonomyCommand,
    ExploreContainmentCommand,
    FindByTaxonomyCommand,
)

__all__ = [
    "UseBranchwaterCommand",
    "ExtractBranchwaterMetadataCommand",
    "ParseContainmentCommand",
    "PlotContainmentCommand",
    "DownloadMetadataCommand",
    "ParseMetadataCommand",
    "CheckMetadataAttributesCommand",
    "CountMetadataCommand",
    "PlotMetadataCountsCommand",
    "SingleSampleCommand",
    "DownloadSraCommand",
    "StatusCommand",
    "ResultsTableCommand",
    "ExtractTargetReadsCommand",
    "DownloadTestGenomeCommand",
    "GenomeSearchCommand",
    "GenomeDownloadCommand",
    "GenomePrepareCommand",
    "EnrichTaxonomyCommand",
    "ExploreContainmentCommand",
    "FindByTaxonomyCommand",
]
