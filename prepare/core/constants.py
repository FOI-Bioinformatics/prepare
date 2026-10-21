"""
Core constants for PREPARE.

This module centralizes all magic numbers and hard-coded values used throughout the application.
"""

# SRA Accession Validation
SRA_VALID_PREFIXES = ("SRR", "ERR", "DRR")
SRA_ACCESSION_PATTERN = r"^[A-Z]{3}[0-9]+$"
GENOME_ACCESSION_PATTERN = r"^GC[AF]_\d{9}\.\d+$"
GENOME_ACCESSION_PREFIXES = ("GCF_", "GCA_")

# Default Thresholds
DEFAULT_CONTAINMENT_THRESHOLD = 0.1
DEFAULT_METADATA_THRESHOLD = 0.0
DEFAULT_SINGLE_SAMPLE_THRESHOLD = 0.1
DEFAULT_STEP_SIZE = 0.1

# Default File Names and Paths
DEFAULT_GENOME_FOLDER = "genomes"
DEFAULT_MATCHES_FOLDER = "matches"
DEFAULT_METADATA_FOLDER = "metadata"
DEFAULT_FASTQ_FOLDER = "fastq"
DEFAULT_PARSED_CONTAINMENT_FILE = "parsed_containment.txt"
DEFAULT_SUMMARY_CONTAINMENT_FILE = "top_containments.txt"
DEFAULT_METADATA_TABLE_FILE = "metadata_table.txt"
DEFAULT_METADATA_COUNTS_FILE = "metadata_counts.txt"
DEFAULT_BRANCHWATER_METADATA_FILE = "branchwater_metadata.txt"
FAILED_ACCESSIONS_FILE = "failed_accessions.txt"

# Threading and Performance Defaults
DEFAULT_NUM_THREADS = 4
DEFAULT_TOP_N = 100
# prefetch's --max-size: a run whose .sra archive is larger is not fetched. The
# prefetch_max_size setting (PREPARE_PREFETCH_MAX_SIZE, config [runtime]) overrides it.
DEFAULT_PREFETCH_MAX_SIZE = "100G"

# Visualization Defaults
DEFAULT_PLOT_COLUMN = "max_containment"
DEFAULT_PLOT_TYPE_CONTAINMENT = "rank"
DEFAULT_PLOT_TYPE_METADATA = "bar"
SUPPORTED_PLOT_FORMATS = ["png", "jpg", "pdf", "svg"]
DEFAULT_FIGURE_SIZE = (12, 8)
DEFAULT_DPI = 300

# Plot Type Choices
CONTAINMENT_PLOT_TYPES = ["rank", "histogram", "box", "violin"]
METADATA_PLOT_TYPES = ["bar", "pie", "radar"]

# File Formats
SUPPORTED_FILE_FORMATS = ["branchwater"]

# Logging Configuration
LOG_LEVELS = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
DEFAULT_LOG_LEVEL = "INFO"

# Security Constraints
# No limit by default; prepare.core.settings' subprocess_timeout setting (--timeout,
# PREPARE_TIMEOUT, config [runtime] timeout) resolves the value prepare.utils.security's
# run_secure actually uses.
DEFAULT_SUBPROCESS_TIMEOUT = 0
# Seconds allowed for a `--version` probe of an external tool: always bounded, regardless of
# the run's own subprocess_timeout setting, so a hung tool cannot stall version detection.
VERSION_PROBE_TIMEOUT = 30
DANGEROUS_ENV_VARS = ["LD_PRELOAD", "LD_LIBRARY_PATH", "DYLD_INSERT_LIBRARIES"]

# Bioinformatics Tool Configuration
ALLOWED_BIOINFORMATICS_TOOLS = {
    "fasterq-dump": {
        "safe_params": {
            "--threads",
            "--progress",
            "-O",
            "--temp",
            "--split-files",
            "--split-3",
            "--skip-technical",
            "--include-technical",
            "--force",
            "--version",
        },
        "description": "NCBI SRA data download tool",
    },
    "prefetch": {
        "safe_params": {
            "-O",
            "--max-size",
            "--progress",
            "--resume",
            "--version",
        },
        "description": "NCBI SRA prefetch (downloads the .sra archive)",
    },
    "pigz": {
        "safe_params": {
            "-p",
            "-f",
            "-k",
            "--version",
        },
        "description": "Parallel gzip",
    },
    "megahit": {
        "safe_params": {
            "-1",
            "-2",
            "-r",
            "-o",
            "--num-cpu-threads",
            "--memory",
            "--min-contig-len",
            "--k-min",
            "--k-max",
            "--k-step",
            "--no-mercy",
            "--bubble-level",
            "--presets",
            "--tmp-dir",
            "--version",
        },
        "description": "Illumina assembly tool",
    },
    "flye": {
        "safe_params": {
            "--nano-raw",
            "--nano-corr",
            "--nano-hq",
            "--pacbio-raw",
            "--pacbio-corr",
            "--pacbio-hifi",
            "--out-dir",
            "--genome-size",
            "--iterations",
            "--meta",
            "--polish-target",
            "--threads",
            "--min-overlap",
            "--keep-haplotypes",
        },
        "description": "Long-read assembly tool",
    },
    "datasets": {
        "safe_params": {
            "download",
            "genome",
            "accession",
            "--include",
            "--filename",
            "--inputfile",
            "--assembly-level",
            "--version",
            "--no-progressbar",
        },
        "description": "NCBI datasets CLI for genome downloads",
    },
    "minimap2": {
        "safe_params": {
            "-a",
            "-x",
            "-t",
            "-o",
            "-d",
            "--secondary",
            "--sam-hit-only",
            "--version",
        },
        "description": "Read-to-reference aligner (targeted read extraction)",
    },
    "samtools": {
        "safe_params": {
            "view",
            "fastq",
            "cat",
            "sort",
            "coverage",
            "-b",
            "-F",
            "-f",
            "-c",
            "-q",
            "-o",
            "-0",
            "-1",
            "-2",
            "-s",
            "-@",
            "--version",
        },
        "description": "SAM/BAM utilities (filter and export mapped reads, reference coverage)",
    },
    "seqkit": {
        "safe_params": {
            "stats",
            "-T",
            "-j",
            "--version",
        },
        "description": "FASTQ/FASTA statistics (fast exact per-file counts for the shared stats cache)",
    },
}

# URL Constants
BRANCHWATER_BASE_URL = "https://branchwater.jgi.doe.gov/"

# Regular Expressions
UNSAFE_SHELL_CHARS = [
    "$",
    "`",
    ";",
    "|",
    "&",
    ">",
    "<",
    "*",
    "?",
    "[",
    "]",
    "(",
    ")",
    "{",
    "}",
]

# File Extensions
FASTQ_EXTENSIONS = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
ASSEMBLY_EXTENSIONS = [".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz"]

# FASTA file patterns accepted as genome inputs (plain and gzipped).
GENOME_FASTA_GLOBS = ("*.fna", "*.fna.gz", "*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz")

# FASTQ file patterns accepted as read inputs (plain and gzipped).
FASTQ_GLOBS = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")

# Most screening entries the registry keeps per target genome; the lowest containments
# beyond this are dropped, so a broad search does not bloat the project journal.
DEFAULT_REGISTRY_MAX_SCREENED = 5000

# Shared data store: one folder holding a single copy of every downloaded
# metagenome, used across organism projects.
STORE_MARKER = "prepare_store.json"
STORE_ENV = "PREPARE_DATA"
CONFIG_DIRNAME = "prepare"
CONFIG_FILENAME = "config.toml"
STORE_LAYOUT = "sra-v1"

# Per-accession dataset lock (prepare.store.locks). A store download or adopt holds its
# accession's lock for as long as the transfer takes, which can be hours, so the holder
# refreshes the lock file's mtime every LOCK_HEARTBEAT_SECONDS and a waiter only reclaims a
# lock whose mtime is older than DATASET_LOCK_STALE_SECONDS, i.e. one whose holder has died.
LOCK_HEARTBEAT_SECONDS = 10.0
DATASET_LOCK_STALE_SECONDS = 600.0

# store_gc keeps a dataset for this long after prepare.store.locks.touch_dataset_use last
# recorded a project starting to use it (a link, a store_link, an adopt), even when that
# project's own usage row has not (yet, or ever) been recorded in the catalogue: a dataset
# just handed to a project must not be removed out from under it by a concurrent gc run.
GC_RECENT_USE_GRACE_SECONDS = 86400.0

# Short-lived locks (prepare.utils.lockfile): the project registry and the store catalogue.
# Their holders write for well under a second, poll quickly while waiting, and refresh the
# lock file every SHORT_LOCK_HEARTBEAT_SECONDS so a slow write is never judged stale. The
# registry's own wait and stale limits live in prepare.data.registry.
SHORT_LOCK_HEARTBEAT_SECONDS = 5.0
SHORT_LOCK_POLL_SECONDS = 0.05
CATALOG_LOCK_WAIT_SECONDS = 60.0
CATALOG_LOCK_STALE_SECONDS = 120.0

# SQLite's own busy handler for every catalogue connection (prepare.store.catalog),
# separate from the O_EXCL lock-file wait above: it covers the brief window between another
# connection taking SQLite's internal lock and catalog_write's own lock file being visible,
# and a store_reindex replay that can hold the database busy longer than sqlite3's implicit
# 5 second default.
CATALOG_BUSY_TIMEOUT_SECONDS = 30.0

# Memory and Resource Limits
MAX_CONCURRENT_DOWNLOADS = 10

# HTML containment explorer: rows written to the results table (highest containment first)
# and the tallest heatmap in pixels, so a large run still opens in a browser.
EXPLORER_TABLE_ROWS = 5000
EXPLORER_HEATMAP_MAX_HEIGHT = 2000

# Plugin System Constants
PLUGIN_REGISTRY_NAME = "prepare_plugins"
