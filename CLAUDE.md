# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PREPARE is a command-line bioinformatics toolkit for analyzing metagenomic datasets based on genome containment. The software processes Branchwater CSV files, downloads SRA metadata from NCBI, and provides advanced visualization capabilities for containment analysis, including intelligent SRA management with resume capability and comprehensive quality profiling.

## Development Commands

### Quick Start
- `make dev-install` - Install with development dependencies
- `make help` - Show all available commands

### Testing and Quality
- `make test` - Run pytest with coverage
- `make lint` - Run flake8 linting
- `make format` - Format code with black
- `make check` - Full quality check (format, lint, type check)
- `make pipeline` - Run full test pipeline via local_test.sh

### Build and Clean
- `make build` - Build distribution packages (uses pyproject.toml)
- `make clean` - Clean build artifacts, cache files, and test data

### Installation
- `make install` - Install package for development
- `make dev-install` - Install with development dependencies
- Legacy: `pip install -r requirements.txt` (core only) or `requirements-all.txt` (core plus every extra)

### Dependencies and extras
- Core runtime: pandas, numpy, matplotlib, biopython, lxml, requests. The CLI must import and build its
  parser with only these (`tests/test_optional_imports.py` blocks every optional module and checks this).
- Extras: `analysis` (scikit-learn, scipy), `interactive` (plotly, jinja2), `maps` (cartopy), `sourmash`;
  `all` installs them all and `dev` depends on `all`.
- Import an optional package at the point of use with `prepare.core.optional.require(module, extra, purpose)`,
  never at module top level and never behind a silent `*_AVAILABLE` flag: a missing package is a
  `ConfigurationError` naming the extra, not a degraded output file. A command that catches broad exceptions
  must re-raise `ConfigurationError` (see `cli/commands/advanced_analysis.py`).
- A test that needs an extra is skipped when it is absent (`pytest.importorskip` or a `requires_*` skip mark),
  so the suite also passes in a core-only install.

## Architecture

PREPARE follows a layered architecture with clear separation of concerns. As of 0.10.0,
`prepare --help` lists 44 commands in seven groups (Containment, Metadata, Genomes, Reads,
Store, Analysis, Environment; the Environment group holds `doctor`, `runs` and `project_report`); the four commands
`sra_profile`/`sra_report` replaced (`sra_stats`, `sra_profile_quality`, `sra_dashboard`,
`sra_compare`) still parse but are hidden from that listing (see "Advanced SRA Commands" below).

### Core Components
- **CLI Layer** (`prepare/cli/`): Command-line interface with modular command architecture using a
  registry pattern. Most commands are one module under `prepare/cli/commands/`; the store and
  status commands are packages instead (`prepare/cli/commands/store/`, one module per command, and
  `prepare/cli/commands/status/`, split into `command.py`, `suggest.py`, `render_text.py`). Every
  command writes to stdout through `BaseCommand.emit`/`emit_raw`/`emit_json` (`cli/base.py`); no
  other module in `prepare/` calls `print(`.
- **Core Logic** (`prepare/core/`): Domain models, validation, exceptions, and constants.
  `core/optional.py` is the one place an optional dependency (scikit-learn, scipy, plotly, jinja2,
  cartopy, sourmash) is imported, at the point of use, raising `ConfigurationError` naming its extra
  when the package is missing. `core/settings.py` resolves every runtime setting (flag,
  `PREPARE_<NAME>`, `config.toml [runtime]`, default); see `docs/configuration.md`.
- **Utilities** (`prepare/utils/`): `tools.py` (the external tool table, `require_tools`),
  `resources.py` (CPUs and memory limit of the job), `progress.py` (progress summaries), `logging.py`,
  `xml.py` (the hardened lxml parser), `http.py` (retrying HTTP sessions), `security.py`
  (`SecureSubprocess.run_secure`), `lockfile.py`, `termination.py`.
- **Data Layer** (`prepare/data/`): File I/O, Branchwater processing, metadata handling, basic SRA
  operations. `data/sra/` is a package (`fastq`, `cleanup`, `accession`, `retry`, `store_handoff`,
  `download`, `space`, `spots`, `run_report`), not a single module; `run_report` builds the
  `download_sra` report CSV and `download_run.json`. `data/registry.py` is the project journal;
  `data/registry_blocks.py` holds the typed dataclass for each block the registry file stores, with
  unknown-key preservation on a round trip. `data/run_log.py` is the per-project run log
  (`<project>/.prepare/runs/`, see "Run log rules" below).
- **Advanced SRA Package** (`prepare/sra/`): Quality profiling and interactive reporting, described
  under "Advanced SRA Commands" below.
- **Processing** (`prepare/processing/`): Containment analysis, statistical processing, counting
  algorithms, and `status_report.py` (builds the `status` command's report from the registry and the
  filesystem, with no dependency on the CLI layer). `project_funnel.py` (the cross-stage funnel of
  `status` and `project_report`), `run_diff.py` (run comparison for `runs`), `project_report.py` (the
  `project_report` dict) and `project_report_markdown.py` (its Markdown rendering) follow the same rule.
- **Plugins** (`prepare/plugins/`): Extensible plugin system for formats and visualizers
- **Visualization** (`prepare/visualization/`): Plotting functionality; `project_report.py` renders
  the `project_report` HTML page (the `interactive` extra). The unused report generator
  `visualization/reporting.py` was removed in 0.9.0.

### Key Design Patterns
- **Command Registry**: All CLI commands inherit from base classes and register themselves via `command_registry`
- **Plugin System**: Uses `PluginRegistry` with automatic discovery for format handlers and visualizers
- **Data Models**: Dataclasses in `core/models.py` define `Containment`, `SRAMetadata`, `GenomeInfo`, and `ContainmentSummary`
- **Dual SRA Architecture**: Basic operations in data layer, advanced features in specialized SRA package

## Project Structure

### Modern Python Packaging
- Uses `pyproject.toml` for modern Python packaging (PEP 518/621)
- Maintains `setup.py` for backward compatibility
- Configuration files: `setup.cfg` (flake8), `pyproject.toml` (build, black, mypy, pytest)

### Documentation Structure
- `README.md` - Main project documentation
- `CONTRIBUTING.md` - Development guidelines
- `docs/ARCHITECTURE.md` - Detailed architecture documentation
- `docs/` - Additional documentation (workflow guides, etc.)
- `CLAUDE.md` - This file for Claude Code guidance

## Important Implementation Notes

### CLI Command Structure
New commands should:
1. Inherit from the appropriate base command class in `cli/base.py`
2. Register themselves in `cli/main.py` via the `register_all_commands()` function
3. Follow the existing pattern for argument parsing and execution
4. **CRITICAL**: Use dashes in CLI arguments (e.g., `--matches-folder`), not underscores

### Advanced SRA Commands
Two analysis commands, both in the "Reads" group, share one statistics path
(`prepare/sra/`: `dataset_stats` finds an accession's mate files and loads the cached statistics
record, `quality` samples per-read figures, `analytics.SRADatasetAnalyzer.profile_dataset_quality`
combines them, `profiles.profile_accession` is the one entry point both commands call):
- `sra_profile` (`cli/commands/sra_profile.py`) - statistics table, one profile JSON per accession,
  registry analysis `"profile"`. Totals and GC come from `store.stats.compute_dataset_stats`; GC is in
  percent (`gc_percent`) everywhere.
- `sra_report` (`cli/commands/sra_report.py`) - one HTML report (`sra_report.html`, figures in
  `sra_report.json`); with `--groups-file` it adds the comparison and statistical tests. Each accession
  is profiled once and the profiles are reused by both sections; registry analysis `"report"`.

`sra_stats`, `sra_profile_quality`, `sra_dashboard` and `sra_compare` (and the dash aliases) were
merged into these in 0.5.0; `cli/commands/renamed.py` keeps the old names as hidden commands that log
the new name and exit 2. `results_table` and `status --export-tsv` read the `"profile"` analysis, or
the old `"sra_stats"`/`"quality"` ones of a registry written before 0.5.0
(`data/registry_blocks.profile_summary`). Shared CLI helpers `read_accessions_file`,
`accessions_from_args` and `resolve_command_store` live in `cli/base.py`.

### Store Commands
`prepare/store/` (package: `resolve`, `layout`, `sidecar`, `catalog`, `link`, `adopt`, `usage`,
`locks`, `stats`, `journal`) and `prepare/cli/commands/store/` implement a shared data store: one
copy of each downloaded SRA accession, reused by every project that links into it. Nine commands,
registered under the "Store" group in `cli/main.py`:
- `store_init` - create a store at `--data-root`, record the project's use of it
- `store_status` - dataset and byte counts, stale projects (`--verbose` also lists every dataset)
- `store_reindex` - rebuild the SQLite catalogue from the sidecar files on disk, replaying the
  append-only journal under `<data-root>/journal/` to restore project and usage records that sidecars
  alone do not carry; if that replay restores no project while the rebuilt catalogue still holds
  datasets, it sets a `rebuilt_without_projects` catalogue flag (cleared once a project is restored,
  by this or a later `store_reindex`); a store folder named like an SRA accession but with no sidecar
  gets one rebuilt from its file names and sizes with state `failed` (pointing at
  `store_verify --rescan --fix-state <ACC>`), a folder whose dataset lock is held is skipped with a
  warning, so is a folder not named like an accession, and an unreadable sidecar still refuses the reindex
- `store_adopt` - fold an existing project `fastq/` folder into the store (`--move` links the project
  folder to the store's copy; `--copy` leaves the project folder as it was, including on a duplicate,
  and still records the copying project as a user of the dataset; an empty accession folder is refused
  either way); the sidecar is written into the staging folder before it is published
- `store_verify` - check a dataset's files against its recorded size, md5 (`--md5`), or NCBI spot
  count (`--spots`, which also falls back to `<data-root>/metadata/` and the calling project's
  `metadata/` folder for an XML `download_metadata` wrote after the fact); `--fix-state` only promotes
  a dataset to `complete` after its files have actually been read through, not on size/md5 alone;
  `--rescan` rebuilds a dataset's recorded file list from what is on disk before checking
- `store_link` / `store_unlink` - add or remove one accession's project symlink; `store_link`, like
  `store_adopt` and `download_sra`, records the verdict `cli/commands/sra_verdicts.store_verdict` gives (the
  sidecar's read count against the project's spot count, merged with the verdict on file), and refuses
  without `--accept-partial` a `partial`, `failed` or sidecar-less copy and one recorded `unverified` but
  short against that spot count
- `store_usage` - which projects and genomes used a given accession, or which datasets are unused
  (`--project`, `--organism`, `--unused`)
- `store_gc` - report (or, with `--yes`, remove) datasets nothing references any more; refuses to run
  when the catalogue holds datasets but no project records at all (run `store_reindex` first), and
  also refuses on a `rebuilt_without_projects` catalogue flag until every project has re-registered
  and `--accept-rebuilt` is given; also lists leftover temp artifacts under `<data-root>/tmp/`, including
  the bare `tmp/<ACC>` staging folder of a killed download once its dataset lock is free; `--json` prints
  the report (or a refusal) as one JSON object

Store discovery rules for agents:
1. **Store required before linking**: `download_sra` only links a dataset and `store_adopt` only
   runs once a store exists, found via `--data-root`, `PREPARE_DATA`, `store.root` in the project
   registry, or `~/.config/prepare/config.toml`. A project that never runs `store_init` (and has
   none of those set) behaves exactly as it did before the store existed; do not assume a store is
   present.
2. **Tests never touch a real store**: every store test uses a `tmp_path` store and monkeypatches
   `PREPARE_DATA` and `XDG_CONFIG_HOME`/`HOME`, so a test run never reads or writes the developer's
   actual store or `~/.config/prepare/config.toml`. Follow this pattern for any new store test.
3. **Explicit-path staging**: `store_adopt` copies a dataset into the store before removing anything
   from the project, so adoption briefly holds up to three copies of one accession (the project's
   original, a staging copy, and the store's copy) and needs at least twice the folder's size free on
   the store's filesystem before it starts (`_has_room_for` in `prepare/store/adopt.py`, through
   `room_shortfall` in `prepare/store/link.py`), only replacing the project folder with a symlink once
   the store copy is verified (and never for `--copy`, which leaves the project folder alone). A
   copy-mode link (`--link-mode copy`) stages the same way in the project: the new copy is made under a
   hidden name after the same free-space check, the earlier copy is moved aside only then, and put back
   whenever the swap does not happen. New code that moves data into or within the store should stage the
   same way rather than renaming in place, so an interruption never leaves a dataset with no complete
   copy anywhere.
4. **Journal replay for usage records**: sidecars record a dataset's own state but nothing about which
   projects used it; that lives only in the SQLite catalogue's `projects`/`usage` tables. `store_reindex`
   rebuilds the catalogue from sidecars and then replays `<data-root>/journal/*.jsonl` (append-only,
   written by every `upsert_project`/`record_usage` call, each line stamped with the time it was
   appended; a project line also carries the hostname that wrote it, a usage line does not) to restore
   those tables. Code that adds a new way to record project or usage data should append to the journal
   the same way, or `store_reindex` will silently lose it.
5. **`rebuilt_without_projects` blocks `store_gc`**: when a `store_reindex` replay restores no project
   at all but the rebuilt catalogue still holds datasets, every dataset would look unused to `store_gc`
   until each project re-registers, so `store_reindex` sets this catalogue flag and `store_gc` refuses
   to run on it. The flag clears itself once a project is restored (by that `store_reindex` or a later
   one); if none ever is, `store_gc --accept-rebuilt` is the only way to unblock it, and only after every
   project using the store has run `store_init` or `store_link` again. Do not treat `--accept-rebuilt`
   as a substitute for that re-registration step.

### Concurrency rules
- Never hold the registry lock across ranking, profiling, read counting, or a store catalogue write;
  do that work against a snapshot loaded without the lock, then record the outcome through
  `registry_update` (`prepare/data/registry_batch.py`) or a `RegistryBatch`, both of which hold the
  lock only for the write itself.
- Write every output file through `prepare/data/file_io.py`'s atomic helpers (`write_text_atomic`,
  `write_bytes_atomic`, `open_atomic`, `write_csv`, `atomic_path`, `copy_file`), never a direct
  `.write_text`, `.write_bytes`, `.to_csv`, `shutil.copy2`/`copyfile`, or `open(path, "w..."/"x...")`;
  `scripts/check_atomic_writes.sh`, part of `make check`, enforces this outside the small allowlist in
  `scripts/atomic_writes_allowlist.txt` (figure writers, `savefig` and plotly `write_html`, are listed
  there: a figure is regenerable and overwritten whole).
- A command that runs external tools on its own worker threads (not the main thread) must pass its
  `Termination.stop` (from `args._termination`, set by `BaseCommand.run`) down to every
  `SecureSubprocess.run_secure` call as `stop=`, the same way the download workers do, or an interrupt
  on the main thread will not stop that thread's tools.
- New lock files reuse `prepare.utils.lockfile.held_lock`/`LockPolicy`, never a hand-written
  `O_EXCL` file; see "Locking" in `docs/ARCHITECTURE.md`.
- Process-level concurrency tests live in `tests/test_concurrency_processes.py` (spawning the real CLI
  with fake tools on `PATH`, no real tool or network); a test that starts a real subprocess carries the
  `multiprocess` pytest marker, which runs by default.

### Run log rules
- A command opts in to the run log by overriding `BaseCommand.records_run(args)` (False by default),
  returning False for an invocation that should leave no record (a dry run, `--no-record`, a read-only
  mode); read flags there with `getattr(args, name, False)`. `cli/main.py` writes the record after the
  command returns; a command never calls `record_run` itself.
- Note figures from `execute` with `prepare.data.run_log.note_run(args, summary=..., detail=...)`, or
  `note_rows(args, rows, *section)` for rows noted one item at a time. Keep the summary small: counts,
  totals and a few settings, never one entry per accession.
- A detail holds rows only: mappings of plain values (no nested mapping inside a row) keyed by the
  accession or by `accession/genome`, under one section per command, for example
  `{"analyses": {"profile": {"SRR1": {...}}}}`. Put no other flat mapping (settings, counts, totals) in a
  detail; `processing/run_diff.py` would read it as a row. Keep field names stable between releases,
  since `runs --diff` and `runs --accession` compare by field name.
- Nothing is written in a folder without a project registry, or when the `run_log` setting is off, and
  a failure to record is a warning that never changes the command's exit code. Do not add a list of runs
  to the registry.

### Settings, errors, logging and tools rules
- **Errors**: a command reports an expected error with `return self.fail(error, "context")`
  (`BaseCommand.fail`, `cli/base.py`), never `self.logger.error(...); return 1`. `fail` logs one line
  with the traceback attached for the log file and returns the exception's exit code: 1 for a
  `PrepareError`, 3 for `ConfigurationError`, 4 for `TransientError` (`NetworkError`,
  `LockTimeoutError`), 130 for an interrupt. Raise the matching class rather than choosing a number;
  exit codes other than 0, 1, 2, 3, 4 and 130 are not used.
- **Per-item logging**: a line logged once per accession, sample or request is logged at
  `prepare.utils.progress.item_level(settings.active().progress_every)` (DEBUG by default), and a long
  loop reports through a `ProgressReporter`; INFO is for summaries and for what a user acts on.
  Warnings and errors about one item keep their level.
- **Settings**: a new tunable value is one entry in `SETTINGS` (`_SPECS` in `core/settings.py`) plus a
  field in `RuntimeSettings`, read with `settings.active().<name>` or, in a command,
  `setting_for(args, name)`; a flag for it uses `default=None` and the spec's `cli_dest`. Never read a
  `PREPARE_*` variable or the config file directly; the store root is the one exception
  (`store/resolve.py`). Add the setting to the table in `docs/configuration.md`.
- **External tools**: a command checks the tools it needs with
  `prepare.utils.tools.require_tools([...])` before any work (skipped for `--dry-run`); a new tool, or
  a new floor, is one entry in `TOOLS` (`utils/tools.py`), which `doctor`, the README's tool table and
  `environment.yml` (checked by `tests/test_environment_pins.py`) must agree with. Tools are run only
  through `SecureSubprocess.run_secure`, whose default timeout is the `subprocess_timeout` setting.

### Plugin Development
- Format plugins inherit from base Plugin class in `plugins/base.py`
- Register with `format_registry` for file format handlers
- Register with `visualizer_registry` for visualization types
- Implement required methods like `validate_header`, `parse_file`, `extract_metadata`

### Data Processing Pipeline
The typical workflow involves:
1. Processing Branchwater CSV files (`use_branchwater`)
2. Extracting/downloading metadata (`extract_branchwater_metadata` or `download_metadata`)
3. Parsing containment data (`parse_containment`) 
4. Visualization and analysis (`plot_containment`, `count_metadata`)
5. Excluding unwanted datasets (`blacklist`)
6. Advanced SRA operations (`download_sra`, `sra_profile`, `sra_report`)
7. Checking project state (`status`)

Each step records its outcome in the project registry (`prepare_registry.json`, see
`prepare/data/registry.py`); `status` reads it back and re-checks the filesystem.

### Code Quality Standards & Current Status

#### Current Quality Status
- **All linting checks passing** - No flake8 violations
- **Code formatting consistent** - Black formatting applied
- **Make check passes** - All quality gates working
- **Test suite stability achieved** - Full test suite passing consistently
- **Runtime warnings eliminated** - Numerical computation warnings resolved
- **API compatibility maintained** - DataFrame deprecation warnings addressed

#### Quality Requirements
Every gate below runs as part of `make check`, in this order; each also runs alone with the
command shown.

- **Formatting**: Black, no reformatting needed; 120 character line length (pyproject.toml).
  `python -m black --check --diff prepare tests`
- **Linting**: flake8 clean. `B902` (bare `except Exception`) and `D101`/`D102`/`D103` (a
  docstring on every public class, method and function) are enforced only on the modules listed
  in `setup.cfg`; every other module has D101-D103 switched off there rather than by a wildcard.
  `python -m flake8 prepare tests`
- **Type hints**: mypy clean. `python -m mypy prepare`
- **No print**: no `print(` or `sys.stdout.write` outside `prepare/cli/base.py`; commands write
  their result through `self.emit`/`emit_raw`/`emit_json`, everything else logs.
  `bash scripts/check_no_print.sh`
- **Plotly CDN**: no reference to the frozen `cdn.plot.ly/plotly-latest` alias; use
  `prepare.utils.html.plotly_cdn_script()` instead.
  `grep -rn "cdn.plot.ly/plotly-latest" prepare --include='*.py'`
- **Complexity**: no function at radon cyclomatic-complexity rank D or worse.
  `python -m radon cc prepare -n D -s`
- **Module size and maintainability**: 800 lines and a maintainability index of 20 or higher per
  module under `prepare/`. `tests/test_module_sizes.py` holds a short, shrink-only
  `KNOWN_EXCEPTIONS` list for modules that were already over a ceiling when the guard was added.
  `python tests/test_module_sizes.py --check`
- **ASCII only**: no non-ASCII byte in `prepare/`, `tests/`, `scripts/`, `Makefile`,
  `setup.cfg` or `pyproject.toml`, aside from the documented exemptions in
  `scripts/check_ascii.sh`: single lines marked `# ascii-ok` with a reason (the package
  author's name, Unicode test fixtures). `bash scripts/check_ascii.sh`
- **Atomic writes**: every file under `prepare/` is written through `write_text_atomic`, `write_csv`,
  `open_atomic` or `atomic_path` (`data/file_io.py`), never `.write_text(`, `.to_csv(<path>` or
  `open(<path>, "w")`, outside `scripts/atomic_writes_allowlist.txt`. `bash scripts/check_atomic_writes.sh`
- **Documented commands**: every command visible in `prepare --help` is mentioned in
  README.md, and README.md invokes no command the CLI registry does not have.
  `python scripts/check_docs_commands.py`
- **Coverage**: minimum 80% for new code, 60% project-wide target; checked by `make test`
  (coverage report), not `make check`.
- **Dependency audit**: `pip-audit --strict` clean; runs in CI (`ci.yml`, `audit.yml`), not in
  `make check`. `pip-audit --strict .`

#### Maintenance Standards
1. **Continuous quality assurance** - Focus on maintaining clean codebase
2. **Add type hints** to remaining modules
3. **Maintain quality standards** with `make check` before commits
4. **Consider pre-commit hooks** for automated quality checks

### Error Handling
- Custom exception hierarchy starts with `PrepareError` in `core/exceptions.py`; each class carries
  its exit code (`ExitCode`: 0, 1, 2, 3, 4, 130), see "Exit-code policy" in `docs/ARCHITECTURE.md`
- Each layer handles errors appropriate to its level
- CLI layer formats errors for user display through `BaseCommand.fail`

## Testing Strategy & Current Status

**Status**: The project maintains **robust test infrastructure** with full test suite stability and comprehensive coverage across core functionality.

### Test Suite Reliability
- **Test Execution**: Full test suite passing consistently with zero failures
- **Runtime Stability**: All numerical computation warnings resolved
- **Test Isolation**: Cross-test contamination issues eliminated
- **Mock Integration**: Proper test doubles configured for external dependencies
- **Edge Case Coverage**: Boundary conditions and error paths comprehensively tested

### Current Coverage Status
- **CLI Commands**: Comprehensive coverage with integration testing
- **Core Processing**: Statistical and diversity analysis modules fully validated
- **Data Layer**: SRA operations and metadata handling thoroughly tested
- **Visualization**: Plot generation and error handling properly covered
- **Advanced Features**: Intelligent download and analytics modules tested

### Testing Commands
- `make test` - Unit tests with coverage reporting
- `make pipeline` - Integration test via `local_test.sh`
- `make check` - Full quality check (format, lint, type check)

`tests/test_performance_regressions.py` holds the project's timing and memory bounds, built on the
synthetic fixtures in `tests/perf_fixtures.py` (a 20,000-dataset registry, a large FASTQ file, a
folder of metadata XML files, and a containment table). Each bound was set to about three times
the time measured on the development machine when the fix it guards landed, so a failure means
either a real regression or a slower machine, not a fixed threshold that was guessed in advance.
The time bounds (not the memory bound) are multiplied by the `PREPARE_PERF_SCALE` environment
variable (default 1, values below 1 are rejected); CI sets it to 4 for its shared runners. The bounds
run by default as part of `make test`; run just these checks with `make test-perf`
(`python -m pytest -m perf`).

### Test Infrastructure Components

#### Core Testing Achievements 
- **Statistical Computing**: Numerical edge cases properly handled
- **API Compatibility**: Future-proof implementation with current libraries
- **Error Propagation**: Comprehensive exception handling across modules
- **Test Isolation**: Module contamination issues resolved

#### Reliability Improvements Implemented
1. **Numerical Stability**:
   - Empty array validation in statistical functions
   - Precision loss detection for nearly-identical data
   - Graceful handling of degenerate cases in PCA and PERMANOVA

2. **Mock Architecture**:
   - Proper cleanup of test doubles to prevent contamination
   - Systematic validation of external API expectations
   - Robust handling of concurrent operations in download tests

3. **Production Bug Prevention**:
   - String replacement edge cases identified and resolved
   - Parameter validation enhanced across SRA operations
   - Function signature consistency verified

#### Test Coverage Distribution
- **Core Functionality**: Comprehensive coverage with edge case validation
- **CLI Interface**: Full integration testing with realistic workflows
- **Data Processing**: Statistical methods and file operations thoroughly tested
- **Visualization**: Plot generation and error handling properly covered
- **Advanced Features**: Concurrent downloads and analytics modules validated

### Testing Guidelines

#### Writing Tests
- **Mock External Dependencies**: Use pytest fixtures for file I/O, network calls
- **Test Error Paths**: Include validation failures, missing files, network errors
- **Edge Cases**: Empty datasets, malformed inputs, boundary conditions
- **Integration**: Test complete workflows with realistic data

#### Test Structure
```python
# Example test pattern for CLI commands
@pytest.fixture
def mock_file_operations():
    with patch('prepare.data.file_io.read_file') as mock:
        yield mock

def test_command_execution(mock_file_operations):
    # Setup mocks
    # Execute command
    # Verify behavior
```

#### Coverage Requirements
- **New Code**: 80% minimum coverage
- **Bug Fixes**: Must include tests for the fixed scenario
- **Critical Paths**: CLI commands, data processing must have high coverage

### Integration Testing

The project includes comprehensive integration testing:
1. `local_test.sh` - Basic end-to-end CLI testing with sample data
2. `tests/test_integration_simple.py` - 12 integration tests covering:
   - SRA metadata workflows (retrieval -> CSV export -> visualization)
   - FASTQ processing pipelines (statistics -> reports -> visualization)
   - Multi-sample comparison workflows
   - Data export to multiple formats
   - Error handling in complete workflows

`tests/test_concurrency_processes.py` runs the CLI as separate processes to check concurrent downloads
into one project and through one store, registry writers during a download, SIGTERM handling, takeover
of a killed lock holder and `store_gc` against a running download. `tests/helpers_processes.py` supplies
the isolated environment (`cli_env`, with `HOME` and the store under `tmp_path`) and fake `fasterq-dump`,
`prefetch` and `pigz` scripts that wait on barrier files, so no real tool or network is needed. These
tests carry the `multiprocess` marker, run by default and are skipped on Windows.

### Quality Gates
Before committing code:
1. `make check` must pass (linting, formatting, type checking)
2. `make test` must pass with coverage improvement
3. `make pipeline` integration test must pass
4. No decrease in overall coverage percentage

### Continuous Integration
Besides the pull-request checks, a nightly GitHub Actions workflow (`.github/workflows/smoke.yml`)
runs `scripts/smoke_chain.sh` against real NCBI/SRA services with conda-installed tools, the only
CI job that touches the network; the same chain runs locally via `make test-network`.

## Development Workflow & Implementation Priorities

### Recommended Development Approach

When working on PREPARE, follow this priority order:

#### 1. Maintenance Tasks (Ongoing Priority)
- **Maintain quality standards**: Run `make format` and ensure `make check` passes
- **No new code without tests**: Continue established testing practices
- **Architecture consistency**: Follow established patterns for new features

#### 2. Current Focus Areas (Optional Enhancement)
- **Remaining visualization modules**: interactive.py, plots.py
- **Processing enhancements**: containment.py optimizations, diversity.py extensions
- **Plugin system expansion**: Additional format handler testing
- **Documentation**: Add testing guides and workflow documentation

#### 3. Code Quality Improvements (Medium Priority)
- **Add type hints**: Complete type annotation for remaining modules
- **Documentation enhancement**: Add docstrings for complex functions
- **Performance optimization**: Large dataset handling and memory efficiency

#### 4. Feature Enhancements (Lower Priority)
- **Only after maintaining current quality levels**
- **Must include comprehensive tests**
- **Should follow existing architectural patterns**

### Implementation Guidelines for Claude

#### When Adding New Features
1. **Check existing patterns**: Look at similar implementations first
2. **Write tests first**: TDD approach preferred for new features  
3. **Mock dependencies**: Don't rely on external systems in tests
4. **Follow plugin patterns**: Use existing registries for extensibility

#### When Fixing Bugs
1. **Write a test that reproduces the bug**: This should fail initially
2. **Fix the bug**: Minimal change approach
3. **Verify the test passes**: And doesn't break existing tests
4. **Check edge cases**: Consider similar scenarios

#### When Refactoring
1. **Ensure full test coverage first**: Don't refactor untested code
2. **Small incremental changes**: One logical change at a time
3. **Maintain backward compatibility**: Especially for CLI interfaces
4. **Update documentation**: If public interfaces change

### File-Specific Guidance

#### Well-Tested Files (Reference Implementations)
- `prepare/cli/commands/*.py` - comprehensive CLI testing patterns
- `prepare/cli/commands/sra_profile.py`, `sra_report.py` - SRA analysis commands
- `prepare/data/file_io.py` - robust file operations
- `prepare/data/branchwater.py` - format handling exemplar
- `prepare/data/metadata.py` - external API integration
- `prepare/data/sra_metadata.py` - XML parsing and NCBI API
- `prepare/data/taxonomy.py` - taxonomy validation and NCBI integration
- `prepare/processing/statistics.py` - statistical analysis
- `prepare/plugins/visualizers/bar.py` - bar chart visualization
- `prepare/sra/reporting.py` - interactive dashboard generation
- `prepare/sra/` - comprehensive test coverage for advanced SRA features

**Test Files**: Extensive test suites in `tests/test_*_extended.py`, `tests/test_integration_simple.py`, and `tests/test_performance_simple.py` provide patterns for comprehensive testing with mocking, fixtures, and edge cases.

#### Remaining Enhancement Opportunities
- `prepare/visualization/plots.py` - Plotting functionality
- `prepare/visualization/interactive.py` - Interactive plotting
- `prepare/processing/containment.py` - Core containment algorithms
- `prepare/processing/diversity.py` - Diversity analysis

### Success Metrics for Development Work

#### Recently Completed
- [x] **Intelligent SRA package implemented** - Complete next-generation SRA capabilities
- [x] **All linting issues resolved** - Clean, consistent codebase
- [x] **CLI commands fully tested** - Coverage raised from an untested baseline to comprehensive
- [x] **Core processing tested** - Coverage raised from an untested baseline to comprehensive
- [x] **Data layer testing completed** - Key modules thoroughly covered
- [x] **Test coverage improvement session** - A large batch of tests added across multiple files
- [x] **Critical modules improved** - sra_reporting, the SRA analysis commands, sra_metadata, bar visualizer, taxonomy
- [x] **Integration test suite created** - End-to-end workflow tests added
- [x] **Performance benchmarks established** - Benchmarked tests added with pytest-benchmark
- [x] **Overall project coverage improved** - Raised substantially from an early baseline
- [x] **Orphan code removed** - Clean architecture maintained

#### Advanced SRA Features Achievements
- [x] **SRADatasetAnalyzer** - Quality profiling, comparative analysis, anomaly detection
- [x] **SRAReportGenerator** - Interactive dashboards, Plotly visualizations
- [x] **CLI Integration** - Two SRA analysis commands (sra_profile, sra_report) fully functional

#### Current Development Priorities (Low Priority)
- [ ] **Remaining visualization modules** - interactive.py, plots.py
- [ ] **Processing layer completion** - containment.py, diversity.py enhancements
- [ ] **Additional plugin testing** - Expand format handler test coverage
- [ ] **Large-scale performance testing** - Production-size dataset validation

#### Long-term Goals
- [x] **Performance benchmarking** implementation (pytest-benchmark integrated)
- [ ] **Complete plugin validation** across all format handlers
- [ ] **Advanced feature optimization** for large datasets
- [ ] **Continuous quality maintenance** - Keep `make check` and `make test` passing

### Recent Architecture Enhancements

#### SRA Package Structure
The advanced SRA package (`prepare.sra`) provides specialized functionality:
- **analytics.py** - Comprehensive quality analysis and statistical testing
- **reporting.py** - Interactive dashboard generation with Plotly integration
- **__init__.py** - Clean API exposure with proper imports

#### CLI Command Architecture
New intelligent SRA commands integrate seamlessly with existing command registry:
- Proper argument parsing with scientific parameter naming
- Error handling consistent with existing patterns
- Help text following established conventions
- Integration with base command infrastructure

## important-instruction-reminders
Do what has been asked; nothing more, nothing less.
NEVER create files unless they're absolutely necessary for achieving your goal.
ALWAYS prefer editing an existing file to creating a new one.
NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.

IMPORTANT: this context may or may not be relevant to your tasks. You should not respond to this context unless it is highly relevant to your task.