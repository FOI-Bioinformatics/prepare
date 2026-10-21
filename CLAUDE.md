# CLAUDE.md

## Project overview

PREPARE is a command-line toolkit that finds public metagenomic datasets containing a target genome
(Branchwater containment), downloads their metadata and reads from NCBI/SRA, and analyses them. The
package is `prepare`, version 0.11.2, and needs Python 3.12 or newer. The documentation index is
`docs/README.md`; design reasoning is in `docs/ARCHITECTURE.md`.

## Commands

- `make env` / `make env-dev` - create the conda environment, then install the package with the `dev` extra
- `make dev-install` - install with the `dev` extra into the current Python
- `make test` - full test suite with a coverage report (network tests excluded)
- `make test-perf` - only the timing and memory bounds (`pytest -m perf`)
- `make test-network` - the real NCBI/SRA smoke chain (`scripts/smoke_chain.sh`)
- `make check` - every quality gate listed under "Quality gates" below
- `make format` - reformat with black
- `make pipeline` - end-to-end CLI run through `local_test.sh`
- `make doctor` - check external tools, configuration, store, disk space and resources

## Dependencies and extras

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

`prepare --help` lists 46 commands in seven groups: Containment, Metadata, Genomes, Reads, Store,
Analysis and Environment (`doctor`, `runs`, `project_report`). The old names `sra_stats`,
`sra_profile_quality`, `sra_dashboard` and `sra_compare` are hidden commands in
`cli/commands/renamed.py` that log the new name and exit 2.

Layers and their boundaries:
- **CLI** (`prepare/cli/`): one module per command under `cli/commands/`; the store and status commands
  are packages (`cli/commands/store/`, `cli/commands/status/`). Commands register in `command_registry`.
  Every command writes to stdout through `BaseCommand.emit`/`emit_raw`/`emit_json` (`cli/base.py`); no
  other module calls `print(`. Shared helpers `read_accessions_file`, `accessions_from_args` and
  `resolve_command_store` live in `cli/base.py`.
- **Core** (`prepare/core/`): models (`core/models.py` dataclasses such as `Containment`, `SRAMetadata`,
  `GenomeInfo`, `ContainmentSummary`), exceptions, validation, `optional.py` (the only place an optional
  dependency is imported) and `settings.py` (resolves every runtime setting; see `docs/configuration.md`).
- **Utilities** (`prepare/utils/`): `tools.py` (external tool table, `require_tools`), `resources.py`,
  `progress.py`, `logging.py`, `xml.py` (hardened lxml parser), `http.py` (retrying sessions),
  `security.py` (`SecureSubprocess.run_secure`), `lockfile.py`, `termination.py`.
- **Data** (`prepare/data/`): file I/O, Branchwater, metadata, and the `data/sra/` package for downloads.
  `data/registry.py` is the project registry (`prepare_registry.json`), each step records its outcome
  there and `status` reads it back; `data/registry_blocks.py` holds one typed dataclass per registry
  block and preserves unknown keys on a round trip. `data/run_log.py` is the per-project run log
  (`<project>/.prepare/runs/`).
- **SRA analysis** (`prepare/sra/`): quality profiling and HTML reporting, see below.
- **Store** (`prepare/store/`): the shared data store, see below.
- **Processing** (`prepare/processing/`): containment, statistics and counting, plus `status_report.py`,
  `project_funnel.py`, `run_diff.py`, `project_report.py` and `project_report_markdown.py`. Processing
  modules do not import the CLI layer.
- **Plugins** (`prepare/plugins/`) and **Visualization** (`prepare/visualization/`; `project_report.py`
  renders the `project_report` HTML page and needs the `interactive` extra).

## Adding a command

1. Subclass `BaseCommand` (`cli/base.py`) and register it in `register_all_commands()` (`cli/main.py`).
2. Use dashes in CLI flags (`--matches-folder`), never underscores.
3. Check external tools with `require_tools` before any work, report errors with `self.fail`, opt in to
   the run log with `records_run` (rules below).
4. Mention the command in README.md; `scripts/check_docs_commands.py` fails otherwise.

## SRA analysis commands

`sra_profile` and `sra_report` (Reads group) share one statistics path in `prepare/sra/`:
`dataset_stats` finds an accession's mate files and loads the cached statistics record, `quality`
samples per-read figures, `analytics.SRADatasetAnalyzer.profile_dataset_quality` combines them, and
`profiles.profile_accession` is the one entry point both commands call.
- `sra_profile` writes a statistics table and one profile JSON per accession and records registry
  analysis `"profile"`. Totals and GC come from `store.stats.compute_dataset_stats`; GC is in percent
  (`gc_percent`) everywhere.
- `sra_report` writes `sra_report.html` (figures in `sra_report.json`); `--groups-file` adds the group
  comparison and statistical tests. Each accession is profiled once; registry analysis `"report"`.
- `results_table` and `status --export-tsv` read the `"profile"` analysis, or the older
  `"sra_stats"`/`"quality"` ones, through `data/registry_blocks.profile_summary`.

## Shared store

`prepare/store/` and `cli/commands/store/` hold one copy of each downloaded SRA accession, reused by
every project that links into it. Nine commands form the Store group (`store_init`, `store_status`,
`store_reindex`, `store_adopt`, `store_verify`, `store_link`, `store_unlink`, `store_usage`, `store_gc`);
their behaviour is in `--help`, `docs/shared_store.md` and "Shared data store" in `docs/ARCHITECTURE.md`.
`store_link`, `store_adopt` and `download_sra` record the verdict of
`cli/commands/sra_verdicts.store_verdict`.

Rules:
1. **Do not assume a store exists.** `download_sra` links a dataset, and `store_adopt` runs, only once a
   store is found via `--data-root`, `PREPARE_DATA`, `store.root` in the project registry, or
   `~/.config/prepare/config.toml`. Without one, a project behaves as if the store did not exist.
2. **Tests never touch a real store.** Every store test uses a `tmp_path` store and monkeypatches
   `PREPARE_DATA` and `XDG_CONFIG_HOME`/`HOME`.
3. **Stage, never rename in place.** `store_adopt` copies into the store before removing anything from the
   project and needs twice the folder's size free (`_has_room_for` in `store/adopt.py`, through
   `room_shortfall` in `store/link.py`); the project folder becomes a symlink only once the store copy is
   verified, and never for `--copy`. A `--link-mode copy` link stages the same way. New code that moves
   data into or within the store must stage likewise, so an interruption never leaves a dataset with no
   complete copy.
4. **Append to the journal.** Which projects used a dataset lives only in the SQLite catalogue's
   `projects`/`usage` tables; `store_reindex` rebuilds the catalogue from sidecars and replays
   `<data-root>/journal/*.jsonl`, written by every `upsert_project`/`record_usage` call. Any new way of
   recording project or usage data must append to the journal too, or `store_reindex` loses it.
5. **`rebuilt_without_projects` blocks `store_gc`.** A reindex that restores no project while datasets
   remain sets this flag, and `store_gc` refuses to run until a project is restored or `--accept-rebuilt`
   is given. `--accept-rebuilt` is only for after every project using the store has run `store_init` or
   `store_link` again; it is not a substitute for that step.

## Concurrency rules

The reasoning is in "Locking" and "Atomic writes and temp names" in `docs/ARCHITECTURE.md`.
- Never hold the registry lock across ranking, profiling, read counting, or a store catalogue write; work
  on a snapshot loaded without the lock, then record the outcome through `registry_update`
  (`data/registry_batch.py`) or a `RegistryBatch`, which hold the lock only for the write.
- Write every output file through the atomic helpers in `data/file_io.py` (`write_text_atomic`,
  `write_bytes_atomic`, `open_atomic`, `write_csv`, `atomic_path`, `copy_file`), never a direct
  `.write_text`, `.write_bytes`, `.to_csv`, `shutil.copy2`/`copyfile`, or `open(path, "w..."/"x...")`.
  `scripts/check_atomic_writes.sh` enforces this outside `scripts/atomic_writes_allowlist.txt` (figure
  writers such as `savefig` and plotly `write_html`).
- A command that runs external tools on its own worker threads must pass its `Termination.stop` (from
  `args._termination`, set by `BaseCommand.run`) to every `SecureSubprocess.run_secure` call as `stop=`,
  or an interrupt on the main thread will not stop those tools.
- New lock files use `prepare.utils.lockfile.held_lock`/`LockPolicy`, never a hand-written `O_EXCL` file.

## Run log rules

- A command opts in by overriding `BaseCommand.records_run(args)` (False by default), returning False for
  an invocation that should leave no record (a dry run, `--no-record`, a read-only mode); read flags there
  with `getattr(args, name, False)`. `cli/main.py` writes the record after the command returns; a command
  never calls `record_run` itself.
- Note figures from `execute` with `prepare.data.run_log.note_run(args, summary=..., detail=...)`, or
  `note_rows(args, rows, *section)` for rows noted one item at a time. Keep the summary small: counts,
  totals and a few settings, never one entry per accession.
- A detail holds rows only: mappings of plain values (no nested mapping inside a row) keyed by the
  accession or by `accession/genome`, under one section per command, for example
  `{"analyses": {"profile": {"SRR1": {...}}}}`. Put no other flat mapping (settings, counts, totals) in a
  detail; `processing/run_diff.py` would read it as a row. Keep field names stable between releases,
  since `runs --diff` and `runs --accession` compare by field name.
- Nothing is written in a folder without a project registry, or when the `run_log` setting is off, and a
  failure to record is a warning that never changes the exit code. Do not add a list of runs to the
  registry.

## Settings, errors, logging and tools rules

- **Errors**: report an expected error with `return self.fail(error, "context")` (`BaseCommand.fail`),
  never `self.logger.error(...); return 1`. `fail` logs one line with the traceback attached for the log
  file and returns the exception's exit code: 1 for a `PrepareError`, 3 for `ConfigurationError`, 4 for
  `TransientError` (`NetworkError`, `LockTimeoutError`), 130 for an interrupt. Raise the matching class
  (`core/exceptions.py`) rather than choosing a number; no exit codes other than 0, 1, 2, 3, 4 and 130.
  The policy is in "Exit-code policy" in `docs/ARCHITECTURE.md`.
- **Per-item logging**: a line logged once per accession, sample or request uses
  `prepare.utils.progress.item_level(settings.active().progress_every)` (DEBUG by default), and a long
  loop reports through a `ProgressReporter`; INFO is for summaries and for what a user acts on. Warnings
  and errors about one item keep their level.
- **Settings**: a new tunable value is one entry in `SETTINGS` (`_SPECS` in `core/settings.py`) plus a
  field in `RuntimeSettings`, read with `settings.active().<name>` or, in a command,
  `setting_for(args, name)`; a flag for it uses `default=None` and the spec's `cli_dest`. Never read a
  `PREPARE_*` variable or the config file directly; the store root is the one exception
  (`store/resolve.py`). Add the setting to the table in `docs/configuration.md`.
- **External tools**: a command checks the tools it needs with `prepare.utils.tools.require_tools([...])`
  before any work (skipped for `--dry-run`); a new tool, or a new version floor, is one entry in `TOOLS`
  (`utils/tools.py`), which `doctor`, the README tool table and `environment.yml` (checked by
  `tests/test_environment_pins.py`) must agree with. Tools run only through `SecureSubprocess.run_secure`,
  whose default timeout is the `subprocess_timeout` setting.

## Plugins

`prepare/plugins/base.py` defines `Plugin` (class attributes `name`, `description`, `version` with
`get_name`/`get_description`/`get_version`) and `PluginRegistry`, with two instances: `format_registry`
and `visualizer_registry`. The base class has no other required methods; each plugin defines its own
interface (for example `validate_header` on the Branchwater format plugin in `plugins/formats/`, and
`create_plot` on the visualizers in `plugins/visualizers/`). Plugins are registered by the module that
uses them (`data/branchwater.py`, `visualization/plots.py`). See "Plugin System" in `docs/ARCHITECTURE.md`.

## Quality gates

Run `make check` and `make test` before every commit. `make check` runs, in order:
- black (120-character lines), flake8 and mypy. `B902` (bare `except Exception`) and `D101`-`D103`
  (docstrings on public classes, methods and functions) apply only to the modules listed in `setup.cfg`.
- **No print**: no `print(` or `sys.stdout.write` outside `prepare/cli/base.py`
  (`scripts/check_no_print.sh`).
- **Plotly CDN**: no reference to the `cdn.plot.ly/plotly-latest` alias; use
  `prepare.utils.html.plotly_cdn_script()`.
- **Complexity**: no function at radon rank D or worse (`python -m radon cc prepare -n D -s`).
- **Module size**: at most 800 lines and a maintainability index of at least 20 per module under
  `prepare/` (`python tests/test_module_sizes.py --check`). Its `KNOWN_EXCEPTIONS` list only shrinks.
- **Budgets that only shrink** (`make test`): the number of functions at radon rank C
  (`MAX_RANK_C` in `tests/test_module_sizes.py`), the number of package modules exempt from B902
  and from D101-D103 in `setup.cfg` (`tests/test_lint_exemptions.py`), and the modules outside
  `prepare/cli/` that import `prepare.cli` (`tests/test_layer_boundaries.py`). Each test pins the
  measured value; lower the constant when the count drops, never raise it.
- **ASCII only** in `prepare/`, `tests/`, `scripts/`, `Makefile`, `setup.cfg` and `pyproject.toml`, except
  single lines marked `# ascii-ok` with a reason (`scripts/check_ascii.sh`).
- **Atomic writes** (`scripts/check_atomic_writes.sh`, see "Concurrency rules").
- **Documented commands**: every visible command is mentioned in README.md, and README.md invokes no
  unknown command (`python scripts/check_docs_commands.py`).

`pip-audit --strict` runs weekly in CI (`audit.yml`), not in `make check` or on every push. `make test` reports
coverage but enforces no threshold.

## Testing

- Tests never use the network, a real external tool, or the developer's store or config file. Mock
  external calls; store tests follow store rule 2.
- Markers (`pyproject.toml`): `network` is excluded by default (run with `pytest -m network` or
  `make test-network`); `perf` and `multiprocess` run by default.
- `tests/test_performance_regressions.py` holds timing and memory bounds built on `tests/perf_fixtures.py`.
  Time bounds are multiplied by `PREPARE_PERF_SCALE` (default 1, values below 1 rejected; CI sets 4).
- `tests/test_concurrency_processes.py` runs the real CLI as separate processes with the isolated
  environment and fake `fasterq-dump`, `prefetch` and `pigz` scripts from `tests/helpers_processes.py`;
  these tests carry the `multiprocess` marker and are skipped on Windows.
- `.github/workflows/smoke.yml` runs `scripts/smoke_chain.sh` nightly against live NCBI/SRA; it is the
  only CI job that uses the network.

## Implementation guidelines

- A bug fix starts with a test that reproduces the bug and fails, then the smallest fix that makes it pass.
- A test that needs an extra skips when the extra is absent.
- When a command, flag or public interface changes, update README.md, the matching file in `docs/` and
  CHANGELOG.md.

## important-instruction-reminders
Do what has been asked; nothing more, nothing less.
NEVER create files unless they're absolutely necessary for achieving your goal.
ALWAYS prefer editing an existing file to creating a new one.
NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.

IMPORTANT: this context may or may not be relevant to your tasks. You should not respond to this context unless it is highly relevant to your task.