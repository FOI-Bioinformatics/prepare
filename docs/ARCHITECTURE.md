# PREPARE Architecture

This document provides an overview of the PREPARE architecture and design principles.

## Design Principles

PREPARE is designed around the following principles:

1. **Separation of concerns**: data access, processing and the command-line interface are separate layers
2. **Extensibility**: format handlers and visualizers are plugins
3. **Consistent commands**: every command inherits from a base class and reports errors and output the same way
4. **Validation at the boundaries**: user input and external data are checked where they enter
5. **Tests for error paths and edge cases**

## Architecture Overview

The application is structured using a layered architecture with the following components:

```
+------------------+
|       CLI        |
+------------------+
         |
+------------------+
|    Core Logic    |
+------------------+
         |
+------------------+     +------------------+
|    Data Layer    |<--->|  Plugin System   |
+------------------+     +------------------+
         |
+------------------+
| External Systems |
+------------------+
```

### Layers

#### CLI Layer
The CLI layer provides the command-line interface for the application. It's responsible for:
- Parsing command-line arguments
- Routing commands to the appropriate handlers
- Formatting output for the user

Most commands are one module each under `prepare/cli/commands/`. Two command groups are large
enough to be their own packages instead: `prepare/cli/commands/store/` (one module per store
command: `init`, `status`, `reindex`, `adopt`, `verify`, `link`, `usage`, `gc`, plus a shared
`_shared.py`) and `prepare/cli/commands/status/` (`command.py` for the `StatusCommand` class,
`suggest.py` for the `--next` suggestions, `render_text.py` for the text-output formatters). Two
commands keep part of their logic in a helper module beside them: `extract_target_reads`
(`read_extraction.py`) assembles through `extraction_assembly.py` (`assemble_samples`, one sample at a
time under its sample lock, with the per-sample outcome in `AssemblyOutcome`), and `download_sra`
(`sra.py`) takes its registry inputs and download verdicts from `sra_verdicts.py` (`registry_inputs`,
`present_verdicts`, `store_verdict`, which `store_link` and `store_adopt` record too) and records
each download outcome in the registry through `sra_record.py` (`result_recorder`, `record_present`).
The two reporting commands of the Environment group are one module each: `runs.py` (`RunsCommand`,
which reads the run log and compares runs through `processing/run_diff.py`) and `project_report.py`
(`ProjectReportCommand`, which builds the report in `processing/project_report.py` and renders it
through `processing/project_report_markdown.py` and `visualization/project_report.py`). Every
command writes to stdout through exactly one channel, `BaseCommand.emit`/`emit_raw`/`emit_json`
(`prepare/cli/base.py`); library modules log or return their output instead of printing, and a
`make check` gate fails on any other `print(` call in `prepare/`.

#### Core Logic Layer
The core logic layer contains the business logic of the application. It includes:
- Data models and domain objects
- Validation logic
- Processing algorithms
- Error handling

#### Data Layer
The data layer manages access to data sources and external systems:
- File I/O operations
- Data format conversions
- Communication with external APIs
- Caching mechanisms

#### Project registry
`prepare/data/registry.py` is the project journal: a `prepare_registry.json` file in the project
root recording screening, selection, exclusion, download, metadata, analysis, extraction, and
assembly outcomes for every accession a project touches. Record functions (`record_screening`,
`record_selection`, `record_exclusion`, `record_download`, `record_metadata`, `record_analysis`,
`record_extraction`, `record_assembly`, `record_genome`) are called by the commands that produce
those outcomes as each one finishes, each inside a `registry_transaction` that loads the file,
applies one record and writes it back under a lock, so a run never reverts another process's edit.
A screening entry carries its own source and query threshold per genome, and the number kept per
genome is capped (`--registry-max-screened`). Presence is never taken from the registry alone:
`status` and the scanning functions (`scan_downloads`, `scan_metadata`, `scan_extractions`,
`scan_assemblies`) re-check the filesystem, so a file removed by hand is reported as missing rather
than done. Writes are atomic (a temporary file renamed into place) and serialized with a lock file
to avoid concurrent corruption. `update_linked` (add or remove one accession from the project's
`store.linked` list) also lives here, next to the other writers.

The registry file carries a `version` field, currently 2. A schema 1 file loads unchanged; any key the
newer schema added is absent and defaults to an empty value. Schema 2 adds two top-level keys: `project`
(`id`, `name`, `path`, `created`, `organisms[]`, `genome_ids[]`), which makes a project self-describing,
and `store` (`root`, `mode`, `linked[]`), which records the shared data store this project uses, if any.
Every path a project registry stores is written relative to the project root when it lies inside that
root, and absolute otherwise, so moving the project directory does not break the registry.

`prepare/data/registry_blocks.py` defines every block the registry file holds (screening,
selection, exclusion, download, metadata, analysis, extraction, assembly, plus the `project` and
`store` blocks above) as a dataclass derived from a common `RegistryBlock`: `from_dict`/`to_dict`
round-trip any key the block does not itself declare, so a field this version of PREPARE does not
know about survives a load-and-save cycle unchanged. `registry.py` reads and writes these typed
blocks rather than raw dictionaries; the on-disk JSON layout is unchanged.

#### Run log
The registry keeps the latest outcome of each step; `prepare/data/run_log.py` keeps one record per
run of a command, so a later run can be compared with an earlier one. The registry itself holds no list
of runs. The layout is under `<project>/.prepare/runs/`, where the project is the folder holding the
registry file:

- `runs.jsonl`: one JSON object per line (`RunRecord`, schema 1), one line per run, in the order the
  runs finished: run ID, command, start and finish times in UTC, seconds, exit code, the argument list
  and parsed arguments with secret values replaced by `***`, PREPARE version, host, process ID, a
  summary and the name of the detail file. Lines are appended under `runs.jsonl.lock` with flush and
  fsync; a line cut short by an interrupted append, or one holding a byte that is not UTF-8 outside a
  JSON string (the file is read with `errors="replace"`), is skipped on reading, with one warning.
- `<run_id>.json`: the detail of one run, written atomically, only when the run noted one. After each
  append, the detail files of that command beyond the last `DETAILS_KEPT_PER_COMMAND` (10) are pruned:
  `runs.jsonl` is first rewritten atomically with `"detail": null` on the pruned lines (each line's own
  JSON is edited, so keys a later schema adds are kept), then the files
  are removed. A pruned run and a run that never noted a detail therefore look the same, and `runs`
  reports both as "not kept".

A command opts in by overriding `BaseCommand.records_run(args)` (False by default) and notes its figures
during `execute` with `run_log.note_run(args, summary=..., detail=...)`, or `note_rows` for rows noted
one item at a time. `cli/main.py` records the run in a `finally` after the command returns, also after
an interrupt (exit 130) or a `SystemExit` (recorded with its own code), when `records_run` is true, the
`run_log` setting is on and a registry file exists (`run_log.project_for`); a failure to record,
including a `records_run` override that raises, is logged as one warning and never changes the exit
code. A detail holds rows only: mappings of plain values keyed by accession or by `accession/genome`,
under one section per kind of result (`download_sra` keeps `failed` and `downloaded`), which is the form
`processing/run_diff.py` compares in `runs --diff` and follows in `runs --accession`. A section is a
mapping whose mapping values are all rows; a flat mapping beside the sections is not a row, and a row
key found in more than one section has its fields prefixed with the section's full path. Readers use
`read_runs` (oldest first; `project_report`'s run section takes the last 10 from it), `resolve_run` (an
ID, a unique prefix, `latest` or `previous`, as `runs` selects a run) and `read_detail`.

#### Optional dependencies
`prepare/core/optional.py` is the single point where an optional package is imported. Core
runtime dependencies (pandas, numpy, matplotlib, biopython, lxml, requests) are always available;
everything else (scikit-learn and scipy behind the `analysis` extra, plotly and jinja2 behind
`interactive`, cartopy behind `maps`, sourmash behind `sourmash`) is imported at the point of use
with `optional.require(module, extra, purpose)`, never at module import time and never behind a
silent availability flag. A missing package raises `ConfigurationError`, naming the extra and the
interpreter to install it into; a command that catches broad exceptions re-raises
`ConfigurationError` rather than reporting a degraded result.

#### Shared data store
`prepare/store` is a package, not a single module, because the store's concerns are independent of
each other and each is small enough to test alone:

- **resolve**: finds the store root, in order, from `--data-root`, `PREPARE_DATA`, `store.root` in
  the project registry, and `[store] data_root` in `~/.config/prepare/config.toml`, and writes the
  `[store]` table of that file (see "Configuration" below; the file reader itself lives in
  `core/settings.py` and is re-exported here).
- **layout**: the store's on-disk shape (`prepare_store.json` marker, `catalog.sqlite`, `locks/`,
  `tmp/`, `sra/<accession>/`) and the paths derived from it.
- **sidecar**: reads and writes `<accession>.json` next to each dataset's files: state, layout,
  compression, per-file size and md5, read and base counts, the NCBI spot and base counts used for the
  completeness verdict, and the cached `stats` block described below.
- **catalog**: the SQLite database (`datasets`, `files`, `projects`, `usage` tables, plus an
  `unused_datasets` view) that lets a query answer "which projects used this accession" or "how many
  bytes belong to this organism" without walking every sidecar. It journals in `DELETE` mode, never
  WAL, since a store is expected to live on a network share and SQLite documents WAL as unsafe there;
  writers are serialised by `catalog.sqlite.lock`.
- **link**: creates and removes the per-accession symlink from a project's `fastq/` folder into the
  store, chooses a relative or absolute target, and detects a dangling link. In copy mode it refreshes an
  earlier copy (a real folder holding `<accession>.json`) by staging: the new copy is made under a hidden
  name, the earlier one is moved aside only then, and it is put back if the swap fails.
- **adopt**: folds an existing per-project `fastq/` folder into the store: copies each accession in
  before removing anything from the project, so a copy always exists somewhere during the operation;
  compares byte content when an accession is already in the store so nothing is duplicated. The sidecar
  is written into the staging folder before it is published, so a published folder always has one.
  `rebuild_missing_sidecar` writes the sidecar `store_reindex` gives a folder that has none: built from
  file names and sizes only, with state `failed` until `store_verify --rescan --fix-state` reads the files.
- **usage**: writes one row per (accession, project, genome, stage) the first and last time each
  combination is used, along with the hostname it ran on, and reports stale projects (whose registry
  is gone, whose `project.id` no longer matches what the catalogue recorded, or that were last seen
  from a different host than the one running the check). On a store shared between machines, a project
  active elsewhere still looks stale from here; `store_gc` leaves a stale project's datasets alone
  unless `--include-stale` is given, and even then a dataset that a stale project still links to is
  kept (the project is listed with " (stale)"). The registry path is stored resolved, so the
  catalogue does not depend on the folder a command was run from.
- **locks**: a per-accession lock file with a heartbeat, so two projects downloading the same accession
  at once cooperate rather than corrupt each other's work; a lock with no heartbeat for 10 minutes
  (`DATASET_LOCK_STALE_SECONDS`) is treated as abandoned and taken over. `--lock-wait` on `download_sra`
  and `store_adopt` bounds how long each waits for another project's lock before giving up.
- **journal**: append-only JSON-lines copies of the catalogue's project and usage records, which
  `store_reindex` replays into a rebuilt catalogue; `store_gc` prunes the lines of the accessions it
  removed and `store_reindex` keeps the latest line per project.
- **sketches**: caches the sourmash signature `sra_sketch` built from a store dataset beside the copy
  (`<ACC>.sig.gz`) and describes it in the sidecar's optional `sketch` block (file snapshot, k-mer
  sizes, totals), with the same file-snapshot invalidation and bounded lock wait as **stats**; the
  hashes stay in the signature file, never in the sidecar.
- **stats**: computes and caches the per-dataset statistics block (streaming exact read and base
  counts, plus a sample for per-read metrics such as GC content), invalidated when the FASTQ file's size
  or modification time changes; used by `sra_profile`, `sra_report` and `sra_validate` (through
  `prepare.sra.dataset_stats.load_dataset_stats` for the first two).

Nothing outside `prepare/store` depends on its internal layout; other layers call its public
functions (`resolve_store_root`, `link_dataset`, `record_usage_safe`, and so on) and otherwise treat a
project with no store configured exactly as it behaved before the store existed.

#### Locking

Every lock file in PREPARE is built on one mechanism, `prepare/utils/lockfile.py`, rather than each
caller managing its own `O_EXCL` file:

- **Holder record**: creating a lock (`O_EXCL`) writes a JSON record of `{pid, host, started, token,
  pidns}` (`token` a random hex string; `pidns` the inode of `/proc/self/ns/pid` on Linux, absent
  elsewhere). The record lets a waiter describe who holds a lock ("pid N on host H since T") and tell a
  dead holder from a live one.
- **Heartbeat**: one daemon thread per process refreshes the mtime of every lock this process holds, at
  an interval set by the lock's policy, so a holder doing slow work is not mistaken for dead.
  `os.register_at_fork` replaces the heartbeat and its thread-local bookkeeping in a forked child, so a
  fork never refreshes or releases the parent's lock.
- **Guarded reclaim**: a waiter that judges a lock stale (its mtime older than the policy's
  `stale_seconds`) or its holder's process dead (same host, `os.kill(pid, 0)` raising
  `ProcessLookupError`) does not unlink it directly. It first creates `<lock>.reclaim` with `O_EXCL`;
  only the process holding that guard re-observes the lock and removes it, and only if the holder record
  and inode still match what was observed and the lock is still stale or dead. This is what keeps two
  waiters judging one lock stale at the same instant from both "winning": only one of them can create
  the guard file.
- **Release**: a lock is removed only by the process whose token matches the one written in the file, so
  a process that took over a lock never removes a later holder's lock on the same path.

Known limitation: a dataset lock's staging folder (a plain project's `.prepare-tmp/<ACCESSION>`, or the
store's own staging path) is named after the accession alone, not the holder's token. A takeover of a
holder that is not dead but merely stalled past the stale window (a suspended process, for example) does
not stop that stalled holder from writing into the same staging folder the new holder is now using, so the
two holders' output can mix there while both run. The mixed folder is never published by the stalled
holder: every destructive or publishing step calls `verify_held` first (the plain-project and store
download publish, the `store_adopt` publish, each registry write and each catalogue commit), and a holder
that finds its lock taken over stops with `LockLost` ("lock lost: ..." for a download) and leaves the
staging folder to the new holder. The per-accession lock file itself is never shared this way; only the
staging path underneath it is.

Four sets of `LockPolicy` limits (`what`, `stale_seconds`, `wait_seconds`, `poll_seconds`,
`heartbeat_seconds`) cover every lock kind:

| Policy | Stale | Wait | Heartbeat |
|---|---|---|---|
| Registry | 120 s | 30 s | 5 s |
| Catalogue | 120 s | 60 s | 5 s |
| Dataset | 600 s | see below | 10 s |
| Run log | 60 s | 10 s | 10 s |

The registry's wait is shorter than its stale window on purpose, and the consequence is visible: a
registry lock orphaned by a holder killed on another host (or in another pid namespace, where it cannot
be judged dead) blocks every registry writer for about 120 s, until the lock goes stale, while each
waiter gives up after its own 30 s and reports the registry as locked. A retry after the 120 s succeeds.
A holder that died on the same host is taken over at once.

The stale, wait and heartbeat values in the table are defaults. Each policy builder reads the runtime
settings (`registry_lock_wait`, `registry_lock_stale`, `catalog_lock_wait`, `dataset_lock_stale`,
`lock_heartbeat`, `lock_wait`) at call time through `settings.setting_or(name, MODULE_CONSTANT)`: a value
the user set by flag, variable or config file wins, otherwise the module constant is used, so a test that
patches the constant still takes effect. A wait that reaches its limit raises `LockWaitTimeout`, a
subclass of both `LockHeld` (so existing handlers keep catching it) and `LockTimeoutError` (exit code 4).

The registry policy covers `<registry>.lock`; the catalogue policy covers
`<store>/catalog.sqlite.lock`; the dataset policy covers four lock files: a store dataset lock
(`<store>/locks/<ACCESSION>.lock`), a plain project's per-accession download lock
(`<fastq>/.locks/<ACCESSION>.lock`), an index build lock (`<index>.lock`), and a per-sample extraction
lock (`<output>/.locks/<ACCESSION>.<GENOME_ID>.lock`). The run-log policy (`RUN_LOG_LOCK_POLICY`,
fixed values, not read from the settings) covers `<project>/.prepare/runs/runs.jsonl.lock`.

The dataset policy's wait varies by caller: a store dataset lock and an index build lock wait without a
time limit (for as long as the holder's heartbeat shows it is alive); a plain project's per-accession
download lock is bounded by `--lock-wait` (0, the default, means no limit); a per-sample extraction lock
never waits (`blocking=False`), so a second process or thread finding a sample already locked reports it
skipped instead of mapping it twice.

Compatibility shims keep every existing caller unchanged: `prepare/store/locks.py` re-exports the
dataset-lock API (`dataset_lock`, `lock_is_held`, `read_holder`, `LockHeld`, `LockWaitStopped`) built on
`held_lock` with the dataset policy, and `prepare/data/registry_io.py`'s `_acquire_lock` is a thin
context manager over `held_lock` with the registry policy.

#### Registry write model

To keep a long command from holding the registry lock across slow work (read counting, profiling,
catalogue writes), the registry is written through one of three patterns, never by loading the file,
doing minutes of work, and writing it back under a single lock:

- **`registry_transaction(path)`**: the low-level context manager, used inside the two helpers below and
  directly by a handful of call sites, that loads the registry under the lock, yields it for one quick
  mutation, and writes it back before releasing the lock.
- **`registry_update(path, mutation)`** (`prepare/data/registry_batch.py`): loads the registry,
  applies `mutation` to it and writes it back, all under one `registry_transaction`, and returns the
  mutation's result. The pattern for a command that works on a read-only snapshot first (ranking,
  parsing, profiling) and then records what it found: `registry_update(path, partial(record_something,
  ...))`. A raising mutation writes nothing.
- **`RegistryBatch`** (`registry_batch(path)`): queues many small mutations (one per accession, for
  example) and applies them in a few transactions instead of one each, flushing on a count or a time
  interval and on `__exit__`; a failing flush is logged rather than raised, since the block's own
  exception, if any, must propagate instead.

`select`, `blacklist`, `download_metadata`, `parse_metadata`, `sra_profile`, `sra_report`,
`sra_validate`, and `status --init`/`--reconcile` all follow this shape: do their own work against a
snapshot loaded without the lock, then record the outcome through `registry_update` or a batch, so the
file they write is re-read under the lock and never overwrites a concurrent `download_sra`'s committed
outcome. `status --reconcile` splits into `registry_reconcile.scan_reconcile` (lists the filesystem once,
no lock, and rehearses the reconcile on a deep copy of the snapshot) and `apply_reconcile` (re-evaluates
the plan's conditions against the registry `registry_update` hands it, reading no files); `reconcile()`
is `apply(scan(...))` for a caller that wants the old one-call behaviour.

#### Download completeness verdicts

A download's verdict compares the reads on disk with NCBI's spot count: `complete` at a ratio of 0.99 or
more, `truncated` below it, `unverified` when no count is known. `prepare/data/sra/spots.py` holds two
rules. `download_sra` (`cli/commands/sra_verdicts.py`, `registry_inputs`), the store hand-off,
`store_link`, `store_adopt` and `status --reconcile` use the lookup, with the project's and then the store's
`metadata/` folder as the XML folders. A plain redownload that comes back `unverified` keeps a recorded
`truncated` verdict, one that comes back `complete` replaces it, and a present download without a known
count is recorded `unverified`. For a store dataset, `sra_verdicts.store_verdict` is the one rule
`download_sra` (a link, or a copy kept but not linked: a result starting with one of
`store_handoff.SETTLED_PREFIXES`), `store_link` and `store_adopt` record: the sidecar's read count judged
against the spot count (`store_spot_count`; the sidecar and the count are read before the registry lock is
taken) and merged with the verdict on file. `store_link` refuses a copy whose sidecar is `unverified` but
short against that count, as it refuses a `partial` one, and rewrites its sidecar as `partial` under the
dataset lock.

- **Spot-count lookup** (`expected_spots(registry, accession, store=None, xml_folders=())`): the first
  positive whole number among the registry's metadata block (`run_total_spots`), the store sidecar's
  `ncbi.spots`, `<ACCESSION>_metadata.xml` in each of the given folders in order, and last the
  `expected_spots` recorded with the previous download verdict. Zero, negative and non-numeric values count
  as unknown. A caller that must not read files (the apply step of a reconcile) passes no store and no XML
  folders.
- **Verdict merge** (`merged_verdict(previous, new, reads_r1, expected)`): when the read count and the spot
  count are both known, the verdict is computed again from them, so it may change in either direction;
  otherwise a missing new verdict keeps the previous one, and a previous `complete` or `truncated` verdict is
  never replaced by `unverified` or by a verdict without a value. A known verdict is never downgraded to
  unknown, since the absence of a count is not evidence about the reads.

In the store, a copy that is `partial` or `failed`, or `unverified` but short against a count known since,
is never linked into a project without `--accept-partial`; such a copy is refetched, and a refetch that is
not `complete` and holds fewer reads than the published copy does not replace it. The sidecar's optional
`refetch` record counts refetches in a row that gained no reads (a read count of 0 is a known count); at
`STORE_PARTIAL_REFETCH_LIMIT` (2, `store_handoff.py`) the copy is no longer fetched without `--force`. A
failure message starting with one of `SETTLED_PREFIXES` (`incomplete:`, `partial in store;`) is a settled
outcome: the retry pass in `data/sra/retry.py` (`_split_not_found`) does not retry it.

#### Assembly identity and staging

`prepare/data/assembly_identity.py` decides whether an assembly folder `<genome_id>_assembly` still
belongs to the reads it would be built from now. A marker file in the folder, `.prepare-assembly.json`
(written with `write_text_atomic`), records `inputs` (the extracted read files as sorted name and size
pairs, the extraction date, the megahit preset, the minimum contig length and the k values), plus the
megahit version, the run parameters and the time written, which are kept for reference and not compared.
The registry's assembly record carries the same `inputs` (`registry_assembly.set_assembly_inputs`).
`assembly_state(out_dir, expected, accept_unmarked)` gives one of four states: `absent` (no folder),
`incomplete` (no `final.contigs.fa`), `stale` (no marker, an unreadable marker, or inputs that differ,
with the differing fields as the reason) or `current`. An unmarked folder of an earlier version is
`current` only when the caller accepts it, which `extract_target_reads` does when the registry record
matches the preset and minimum contig length and is not older than the extraction; it then writes a marker
whose megahit version is empty.

megahit never writes into the final folder. It writes into a hidden staging folder beside it (a
`unique_temp_path` name, `.<name>.<host>.<pid>.<token>.tmp`), the marker is written there, and
`publish_assembly` renames the existing folder aside, renames the staging folder into place, and removes
the aside copy; if the second rename fails, the aside copy is renamed back. A failed run removes its
staging folder, so the previous assembly stays in place. `sweep_staging`, called by `extract_target_reads`
under the sample's extraction lock before it decides on a sample, removes staging and aside folders a
killed run left behind. When the final folder is missing and both an aside copy and a staging folder are
there (a run stopped between the two renames), it first renames a staging folder holding both
`final.contigs.fa` and the marker into place, or else the newest aside copy. An aside copy with no
staging folder beside it is left over from a finished publish, so a final folder missing then was
removed afterwards and is not brought back. Staging and aside names are hidden, so `scan_assemblies` and
`visible_files` never list them.

#### Atomic writes and temp names

`prepare/data/file_io.py` provides `atomic_path`, `write_text_atomic`, `write_bytes_atomic`,
`open_atomic`, and `write_csv`, every one of them writing into a uniquely named temporary file next to
the target (`.{name}.{host}.{pid}.{token}.tmp`, dot-prefixed so it is invisible to every folder listing
built on `visible_files`) and publishing it with `os.replace`, so a reader never sees a partially
written file and a crash leaves at most an orphaned, clearly-named temp file. The registry and store
sidecar writes also call `fsync` before replacing; other outputs do not. `scripts/check_atomic_writes.sh`,
run as part of `make check`, fails on a direct `.write_text`, `.write_bytes`, `.to_csv`, or `open(path,
"w...")` call outside the allowlist in `scripts/atomic_writes_allowlist.txt`. The allowlist holds
`file_io.py` itself and the figure writers (`savefig` and plotly `write_html`, whose output is
regenerable and overwritten whole): `visualization/interactive.py` and `plots.py`, and the bar, heatmap
and map visualizers.

#### Termination

`prepare/utils/termination.py` installs a signal handler for `SIGINT`, `SIGTERM`, and `SIGHUP` (and
`SIGBREAK` in place of the latter two on Windows) for the length of one command, through
`graceful_termination()`. The first signal sets a `Termination.stop` event and raises
`KeyboardInterrupt`, so `finally` blocks and `with` exits (a registry batch flush, a tool shutdown) run
normally; a second signal is logged and counted, not raised, so it cannot cut that flush short; after
`ABANDON_AFTER_SIGNALS` (3) signals in total the process terminates any running child and exits at once
with status 130 (`EXIT_INTERRUPTED`). `BaseCommand.run` (`prepare/cli/base.py`) opens this context
around every command's `execute()`, unless a command sets `graceful_shutdown = False`, puts the
`Termination` on `args._termination`, and on a caught `KeyboardInterrupt` logs, calls
`SecureSubprocess.terminate_children()`, and returns 130.

A per-run stop token (an `Event`, not the process-wide `STOP`) follows a download from `download_sra`
through `_project_download`/`_store_download`, `download_accession`, and `compress_fastq` down to
`SecureSubprocess.run_secure`, so interrupting one run's tools (`terminate_children(stop=token)`) does
not affect a concurrent run's tools in the same process. `STOP` and `SecureSubprocess._stopping` remain
the process-wide emergency stop, set only by the final abandon path and never cleared at the start of an
ordinary run, so an in-process caller that starts a further run after an interrupted one must clear it
itself if it wants that further run's own tools to start.

### Key Components

#### Core Components

- **Models**: Data classes representing domain objects like `Containment` and `SRAMetadata`
- **Validation**: Functions for validating input data and configurations
- **Exceptions**: Custom exception hierarchy for clear error reporting

#### Data Access Components

- **file_io**: Abstract file operations for reading/writing various file formats
- **branchwater**: Functionality for working with Branchwater data
- **metadata**: Functions for downloading and processing metadata
- **sra** (`prepare/data/sra/`): A package, not a single module, for downloading and working with
  SRA data: `fastq` (finding and reading FASTQ files, verifying a download), `cleanup` (transient
  folder handling), `accession` (running `prefetch`/`fasterq-dump` for one accession), `retry`
  (parallel download with retries), `store_handoff` (linking a download into the shared store),
  `spots` (the expected spot count and the verdict merge rule, described in "Download completeness
  verdicts" below), `sampling` (`sample_records`, the index-based FASTQ sampler behind the quality
  profile and the store statistics record, and the Phred histogram helpers), and `download` (the CLI-facing entry point); the package's `__init__.py` re-exports
  the public functions other layers import, so `from prepare.data.sra import download_accession` still
  works
- **accessions** (`prepare/data/accessions.py`): `require_accession` and `read_accessions_file`, the
  one check of an SRA run accession before it becomes a folder name, and the one reader of an
  accessions file (blank lines and `#` comments skipped, repeats dropped, each line checked)
- **registry**: The project journal (`prepare_registry.json`); records dataset state and
  re-checks presence against the filesystem
- **registry_io**: loading, the lock, atomic saves and `registry_transaction`; split out of
  `registry.py` and re-exported from it
- **registry_scan**: the disk scanners (`scan_downloads`, `scan_metadata`, `scan_extractions`,
  `scan_assemblies`) and `bootstrap_from_disk`; split out of `registry.py` and re-exported from it.
  Both modules look `registry` up at call time, since `registry.py` imports them
- **registry_blocks**: Typed dataclasses for every block the registry file holds, described in
  "Project registry" above
- **registry_timing**: `started`/`seconds` for the download, extraction and assembly blocks
  (`set_download_timing` and its two siblings, `Stopwatch`, and `timing_summary` for `status`), kept
  out of `registry.py`, which is at its size ceiling
- **sra/space** (`prepare/data/sra/space.py`): the free-space guard of `download_sra`. `SpaceGuard`
  groups the output, temporary and `.sra` cache folders by filesystem and reserves, per accession,
  10 times its registry run size for output (`FASTQ_EXPANSION` plus `GZIP_EXPANSION`) and 8 times
  for temporary files, or
  `min_free_gb` when the size is unknown. A reservation that does not fit waits for running downloads
  to release theirs; one that could not fit even then fails that accession alone
  (`insufficient-space: ...`). Only a tool's own out-of-space error aborts a download pass
- **assembly**: the megahit step of `extract_target_reads` (`assemble_extracted_reads`,
  `_megahit_args`, `summarise_contigs`), split out of `read_extraction.py` and re-exported from it
- **read_extraction_export**: the minimap2 and samtools argument builders, the export of mapped
  reads to FASTQ and the per-sample reference coverage table; split out of `read_extraction.py` and
  re-exported from it
- **metadata_fetch**: the Entrez efetch calls of `download_metadata` (single and batched, pacing and
  credentials, the split of an efetch response into one package per run); split out of `metadata.py`
  and re-exported from it
- **ncbi_pacing**: `RequestPacer`, the thread-safe request pacing shared by the NCBI taxonomy client,
  `SRAMetadataClient` and the Entrez metadata download: 3 requests per second without an API key and
  10 with one
- **assembly_identity**: the assembly marker, staging and publish (`AssemblyInputs`, `assembly_state`,
  `publish_assembly`, `sweep_staging`), described in "Assembly identity and staging" below
- **registry_assembly**: the `inputs` entry of the registry's assembly record (`set_assembly_inputs`)
  and the two checks that decide whether an unmarked assembly of an earlier version is accepted
  (`legacy_assembly_current`, `assembly_predates_extraction`), kept out of `registry.py`, which is at
  its size ceiling
- **metadata_fields**: the mapping of NCBI metadata columns to registry metadata fields
  (`metadata_fields`, used by the metadata commands in `cli/commands/metadata.py`), and
  `fill_metadata_from_xml`, which fills a metadata block recorded without a spot count from
  `<ACCESSION>_metadata.xml` when that file records a positive count; used by `status --init` and
  `status --reconcile`
- **registry_reconcile**: `scan_reconcile`/`apply_reconcile` for `status --reconcile` (see "Registry
  write model" below); the report, `StoreReconcileReport`, adds `store_unavailable`, `metadata_filled`,
  `verdicts_rechecked` and `assemblies_dropped` (assembly records dated before their extraction record,
  which the apply step removes with `clear_assembly`) to the fields of `ReconcileReport`
- **store** (`prepare/store/`): The shared data store package, described in "Shared data store"
  below; a project that never runs `store_init` never touches it
- **run_log**: the per-project run log (`RunRecord`, `note_run`, `note_rows`, `record_run`, `read_runs`,
  `read_detail`, `resolve_run`, `project_for`), described in "Run log" above
- **sra/run_report** (`prepare/data/sra/run_report.py`): the per-run summary of `download_sra`:
  `failure_reason` (one of `network`, `not-found`, `disk-full`, `insufficient-space`, `locked`,
  `interrupted`, `unknown`), `RunOutcomes` (the last outcome and the number of attempts that started a
  download, per accession), the rows of the `--report-file` CSV and the `download_run.json` document,
  and `run_log_entries`, the run-log summary and detail taken from that document. Not re-exported from
  `prepare/data/sra/__init__.py`

#### Processing Components

- **containment**: Algorithms for analyzing containment data
- **counts**: Functions for counting and summarizing metadata
- **statistics**: Statistical analysis utilities, including `compare_group_means` (t-test for two
  groups, one-way ANOVA for more; used by `sra_report --groups-file`)
- **status_report**: Builds the `status` command's report (present/missing reconciliation, stage
  filtering, store link status) from the registry and the filesystem; kept in `processing/` rather
  than `cli/` so it has no dependency on the CLI layer, and returns data that `cli/commands/status/`
  formats for text or JSON output
- **doctor_report**: the checks of `prepare doctor` (Python, tools, config file and settings, store,
  free space, registry, resources, optionally the network), each a `Check(name, status, detail, data)`;
  `cli/commands/doctor.py` renders them as text or JSON and exits with 3 when one failed
- **project_funnel**: `funnel(registry, members)`, the number of datasets screened, selected,
  downloaded, analysed, extracted and assembled, with the downloaded bytes and time (the time of failed
  downloads apart, as `failed_seconds`) and, for
  extractions and assemblies, the accession and genome pairs that reached the stage with their time and
  assembled bases; the counts are those of `stage_members`, so they match the stages of `status`. Used
  by `status` (the funnel line and the `funnel` key of `--json`) and `project_report`
- **run_diff**: pure functions over run-log records: `diff_summaries` (every summary key of two runs,
  with the numeric difference), `detail_rows` and `diff_details` (rows added, removed and changed between
  two details), `detail_kept` and `accession_history` (the runs whose kept detail holds a row for one
  accession). Used by `runs`
- **project_report**: `build_project_report(registry, max_rows, include_environment, runs_limit)`, the
  report of `project_report` as one dict with ten sections (project, funnel, genomes, extractions,
  downloads, failures, timing, environment, outputs, runs). It reuses `stage_members`, `funnel`,
  `download_verdicts`, `genome_counts`, `timing_summary`, the rows of `results_table`, `failure_reason`
  and the `doctor` checks without network access; it imports nothing from the CLI layer, opens no FASTQ
  file and writes nothing
- **project_report_markdown**: `render_markdown(report)`, and the block layout (text lines and tables)
  that the HTML renderer shares, so both formats show the same content

#### Visualization Components

- **plots**: Functions for generating various types of plots
- **project_report** (`visualization/project_report.py`): `render_html(report)`, the HTML page of
  `project_report`, one `<section>` per report section, with a funnel figure and a per-genome figure;
  plotly and jinja2 are imported through `optional.require` (the `interactive` extra), and plotly.js is
  embedded in the page

#### SRA analysis package

`prepare/sra/` holds the analysis behind `sra_profile` and `sra_report`; the data layer's
`prepare/data/sra/` holds downloading. The modules are:

- **dataset_stats**: `accession_fastq_files` finds the mate files of an accession and
  `load_dataset_stats` loads its statistics record (exact read and base totals and GC), which is cached
  in the store sidecar when the accession is a store link
- **quality**: `SequenceQualityAnalyzer` draws a sample of reads from every mate file, through
  `data/sra/sampling.py`, and returns per-read length, GC, base quality, N content, complexity,
  duplication and adapter figures; base qualities are kept as a Phred histogram
- **analytics**: `SRADatasetAnalyzer` combines the statistics record with the sampled figures
  (`profile_dataset_quality`) and holds the group comparison, anomaly detection and recommendations; the
  result types are `QualityProfile`, `ComparativeAnalysis`, `AnomalyReport` and
  `ProcessingRecommendations`
- **profiles**: `profile_accession`, the one entry point both commands call, plus the profile JSON
  files (`write_profile_json`, `load_quality_profiles`) and `statistics_row`, the row of the
  `sra_profile` table
- **reporting**: `SRAReportGenerator` writes the HTML report; plotly and jinja2 come from the
  `interactive` extra and are checked before any profiling starts

`sra_profile` and `sra_report` share this path, so a dataset is profiled by one set of rules whichever
command runs it.

Other modules, one line each:

- `data/`: `branchwater_search` (API client and cache), `screening_table`, `genome_download`, `gtdb`,
  `genome_taxonomy`, `sra_metadata` (NCBI metadata client), `extraction_locks`, `extraction_staging`
- `processing/`: `diversity`, `results` (the `results_table` rows), `doctor_report`, `selection`
- `visualization/`: `explorer` (`explore_containment`), `interactive`
- `utils/`: `browser`, `html` (including `plotly_js_script`), `sizes` (`tree_bytes`, the size of a
  folder's regular files without following symlinks, used by `store_gc`, `store_link` and cleanup)
- `cli/commands/`: `branchwater_search`, `explore`, `results`, `samples`, `sra_enhanced`,
  `advanced_analysis`, `test_data`
- `plugins/sourmash_plugin.py`: see "Plugin System"

#### Utility Components

- **logging** (`utils/logging.py`): `setup_logging` with a console handler and an optional file
  handler; see "Logging policy" below
- **progress** (`utils/progress.py`): `ProgressReporter`, the thread-safe summary line of a long loop
  (every `progress_every` items and at least every 5 minutes, then a closing line), and `item_level`,
  the level for per-item lines (DEBUG, or INFO when summaries are turned off)
- **tools** (`utils/tools.py`): the one table of external tools, `TOOLS` (tool, conda package, oldest
  supported version, commands that use it, whether it is optional); `probe_tool` finds a tool on
  `PATH` and reads its version with a fixed 30 s timeout; `require_tools(names)` raises one
  `ConfigurationError` (exit 3) listing every tool that is missing or below its floor, with the
  `conda install` command for each. Commands call it before any work; `doctor` reads the same table
- **resources** (`utils/resources.py`): `available_cpus` (affinity mask, then `SLURM_CPUS_PER_TASK`,
  then the CPU count), `memory_limit_bytes` (cgroup v2, cgroup v1, then `SLURM_MEM_PER_NODE` or
  `SLURM_MEM_PER_CPU`; None on macOS or without a limit) and `parse_memory` for `--assembly-memory`.
  The default download worker count and the megahit `--memory` value come from here
- **xml** (`utils/xml.py`): `SAFE_PARSER` and `parse_xml_file`, the lxml parser every lxml call site
  uses, with entity resolution, network access and external DTD loading refused. The standard
  library `ElementTree` call sites use expat, which does not resolve external entities (so XXE does
  not apply) but is not hardened against entity expansion; they parse only small responses fetched
  directly from NCBI
- **http** (`utils/http.py`): `retrying_session`, a `requests.Session` that retries GET requests on
  connection errors and HTTP 429 and 5xx with backoff; used by the NCBI taxonomy, GTDB and
  `SRAMetadataClient` clients. A request that still fails with one of those raises `NetworkError`
- **security** (`utils/security.py`): `SecureSubprocess.run_secure`, the one way an external tool is
  started, with an argument allow-list, child tracking for termination, and the timeout taken from the
  `subprocess_timeout` setting (0, the default, is no limit). On POSIX each tool starts in a session of
  its own, and a timeout or a termination signals its whole process group, so a process the tool
  started (megahit's `megahit_core`) stops with it

## Plugin System

Format handlers and visualizers are plugins, using a registry pattern (`prepare/plugins/base.py`):

1. `Plugin` is an abstract base with the class attributes `name`, `description` and `version`;
   `FormatPlugin` (abstract `validate_header` and `parse_file`) and `VisualizerPlugin` (abstract
   `create_plot`) define the two plugin interfaces
2. `PluginRegistry` has two instances: `format_registry` for file formats and `visualizer_registry` for
   visualization types. `register` refuses an abstract class and a different class under a name already
   taken; registering the same class again does nothing
3. The built-in plugins are registered by `prepare.plugins.load_builtin_plugins()`, which the registries
   call on their first lookup, so importing `prepare.plugins` (as the sourmash entry point does) does
   not import pandas or matplotlib
4. Third-party plugins are found through the `prepare_plugins` entry-point group
   (`discover_entry_points`); one that fails to load is logged at WARNING and skipped

The one format plugin, `BranchWaterFormatPlugin`, lives in `prepare/plugins/formats/`, and the
visualizers (`bar`, `heatmap`, `map`) in `prepare/plugins/visualizers/`. `prepare/plugins/sourmash_plugin.py`
is separate: it exposes the parse, plot, diversity and taxonomy commands as `sourmash scripts`
subcommands.

## Error Handling

Errors follow these rules:

1. **Exception hierarchy**: custom exceptions start at `PrepareError`, with one class per error category
2. **Layered handling**: each layer handles the errors appropriate to its level
3. **Numerical edge cases**: statistical functions check for empty or degenerate input and report it
4. **Input validation**: user input and external data are validated before use
5. **Display**: the CLI layer catches errors and formats them for the user

### Exit-code policy

Each exception class carries the process exit code it gives (`PrepareError.exit_code`,
`core/exceptions.py`), and `exit_code_for(error)` maps an exception to a code:

| Code | `ExitCode` | Raised as |
|---|---|---|
| 0 | `OK` | |
| 1 | `FAILURE` | `PrepareError` and every subclass not below; any other exception |
| 2 | `USAGE` | `ValidationError` (invalid input, such as a malformed accession); argparse usage errors; a renamed command (`cli/commands/renamed.py`) |
| 3 | `CONFIGURATION` | `ConfigurationError` (see below) |
| 4 | `TRANSIENT` | `TransientError`: `LockTimeoutError`, `NetworkError` (see below) |
| 130 | `INTERRUPTED` | `KeyboardInterrupt`, which the first `SIGINT`, `SIGTERM` or `SIGHUP` becomes |

`ConfigurationError` covers a missing optional package, a missing or too old external tool
(`require_tools`), a malformed config file or setting, no NCBI email address, and a log file that
cannot be opened. `LockTimeoutError` is a lock wait that reached its limit; `NetworkError` is a
connection error, a timeout, or HTTP 429 or 5xx after the retries of `utils/http.py`.

A command reports an expected error with `return self.fail(error, "context")` (`BaseCommand.fail`,
`cli/base.py`), which logs one line with the exception attached and returns `exit_code_for(error)`.
`main()` does the same for an exception a command lets through. No code 5 or higher is used.

Two commands decide their code from a set of per-item outcomes: `download_sra` returns 4 only when
every accession that failed was classified `network` by `classify_download_error`, and 1 otherwise
(a truncated download is a failure unless `--accept-partial` is given, as an incomplete store copy
already was).
`download_metadata` logs each accession NCBI did not return, records what it fetched and exits 1
when any accession is missing, or 4 when every failure was a network one and nothing was fetched
(`NetworkError`); a rerun fetches only the missing ones. `parse_metadata` likewise writes the table
from the files it could parse and exits 1 when any file could not be parsed.

### Logging policy

- A command's result goes to stdout through `BaseCommand.emit`/`emit_raw`/`emit_json`; everything else
  is logged to stderr. `make check` fails on any other `print(`.
- `setup_logging` (`utils/logging.py`) installs a console handler and, with `--log-file`, a file
  handler that appends, records host and process ID on every line
  (`%(asctime)s %(hostname)s[%(process)d] %(levelname)s %(name)s: %(message)s`) and receives INFO and
  above even when the console is at WARNING. Only handlers it installed itself are replaced on a
  second call.
- Tracebacks go to the log file only. `ConsoleFormatter` hides the exception text on the console
  unless the level is DEBUG, so an expected failure is one line there.
- Per-item lines (one per accession or request) are logged at `progress.item_level(every)`, which is
  DEBUG; a `ProgressReporter` logs one summary line every `progress_every` items and at least every
  5 minutes at INFO. `--progress-every 0` turns the summaries off and puts the per-item lines back at
  INFO. Warnings and errors about one item keep their level. `data/sra/retry.py`,
  `data/metadata.py`, `data/sra/accession.py`, `data/sra/cleanup.py` and `store/link.py` follow this
  (`progress.active_item_level()` reads the setting). `extract_target_reads` reports from its
  `on_result` callback, and a `progress.DemoteInfo` filter moves the per-sample INFO lines of
  `data/read_extraction.py` to the item level for the length of the run.
- At DEBUG `main()` logs the version, the command line (with the `--api-key` value hidden), host,
  process ID, the SLURM job and array task IDs when set, and every runtime setting with its source.

## Configuration

Two modules resolve configuration, with different orders of precedence (see
[configuration.md](configuration.md) for the user-facing description):

- **Runtime settings** (`core/settings.py`): `SETTINGS` holds one `SettingSpec` per value a run can be
  tuned with (name, parser, default, environment variable, description, flag `dest`, config key). Each
  is resolved from the flag (when given), then `PREPARE_<NAME>` (and any alias such as
  `NCBI_API_KEY`), then the `[runtime]` table of `$XDG_CONFIG_HOME/prepare/config.toml`, then the
  default; `Resolved.source` records which. `main()` calls `settings.activate(args)` once, before
  logging is set up; code reads `settings.active()`, or `setting_for(args, name)` inside a command so a
  flag in a hand-built namespace is honoured. `active()` without `activate()` (library use) resolves
  from the environment and file alone. A value that does not parse, a malformed file, or two lock
  limits that contradict each other raise `ConfigurationError`; `main()` then stops every command
  except `doctor`, which runs on the defaults (`activate_defaults`) and reports the error as a failed
  check. Unknown `[runtime]` keys are kept in `RuntimeSettings.warnings` and logged once logging is up.
- **Store root** (`store/resolve.py`): `resolve_store_root` takes the first of `--data-root`,
  `PREPARE_DATA`, the registry's `store.root` and `[store] data_root` in the same file, and requires
  the store marker. The registry comes before the per-user file because a project that joined a store
  records it, and every command in that project must find the same store; a runtime setting has no
  project-level record. `store_init --set-default` writes the `[store]` table.

`config_path` and `read_config` live in `core/settings.py` and are re-exported from `store/resolve.py`.
