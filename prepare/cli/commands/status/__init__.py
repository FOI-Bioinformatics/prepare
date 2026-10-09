"""
Local inventory / status CLI command.

Reports what PREPARE has already downloaded locally (SRA reads, NCBI metadata,
genome assemblies) so a user can see what is available without re-downloading,
and reports where every accession sits in the registry (screened, selected,
excluded, downloaded, analysed, extracted, assembled). When no registry file
exists yet, the report is reconstructed in memory from what is on disk.

`command` holds `StatusCommand`; the report itself is built by
`prepare.processing.status_report`, the suggested next steps by `suggest` and the text
output by `render_text`.
"""

from prepare.cli.commands.status.command import StatusCommand

__all__ = ["StatusCommand"]
