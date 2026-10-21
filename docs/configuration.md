# Configuration

PREPARE reads its run-time settings from command-line flags, environment variables and an optional
user config file. This page lists every setting, where each one can be given, and the order in which
the sources are consulted. The settings are defined in one table, `SETTINGS` in
`prepare/core/settings.py`; the location of the shared data store is resolved separately, by
`prepare/store/resolve.py`.

## The config file

The file is `$XDG_CONFIG_HOME/prepare/config.toml`, or `~/.config/prepare/config.toml` when
`XDG_CONFIG_HOME` is not set. It is optional. It has two tables:

```toml
[runtime]
timeout = 0                    # seconds before an external tool is stopped; 0 means no limit
min_free_gb = 50
progress_every = 100
log_host = true
ncbi_email = "you@example.org"

[store]
data_root = "/proj/group/prepare-store"
```

`store_init --set-default` writes the `[store]` table and keeps every other table in the file.

A file that is not valid TOML, a `[runtime]` that is not a table, or a value that does not parse
stops every command with exit code 3 and a message naming the file and the key. `prepare doctor` is
the exception: it still runs and reports the problem as a failed check. A key in `[runtime]` that
names no setting is reported as a warning and ignored. Other tables are ignored.

Values may be written as TOML numbers, booleans or strings (`timeout = 600` and `timeout = "600"` are
the same). An array or a table is not a valid value.

## Runtime settings

Each setting can be given in the `[runtime]` table under its config key, in its environment variable,
and, for most, with a flag. An empty environment variable counts as not set.

| Config key | Type | Default | Environment variable | Flag |
|---|---|---|---|---|
| `timeout` | seconds, 0 or more | 0 (no limit) | `PREPARE_TIMEOUT` | `--timeout` |
| `max_workers_cap` | whole number, 1 or more | 4 | `PREPARE_MAX_WORKERS_CAP` | none |
| `lock_wait` | seconds, 0 or more | 0 (no limit) | `PREPARE_LOCK_WAIT` | `--lock-wait` |
| `registry_lock_wait` | seconds, more than 0 | 30 | `PREPARE_REGISTRY_LOCK_WAIT` | none |
| `registry_lock_stale` | seconds, more than 0 | 120 | `PREPARE_REGISTRY_LOCK_STALE` | none |
| `dataset_lock_stale` | seconds, more than 0 | 600 | `PREPARE_DATASET_LOCK_STALE` | none |
| `lock_heartbeat` | seconds, more than 0 | 10 | `PREPARE_LOCK_HEARTBEAT` | none |
| `catalog_lock_wait` | seconds, more than 0 | 60 | `PREPARE_CATALOG_LOCK_WAIT` | none |
| `ncbi_email` | email address | none | `PREPARE_NCBI_EMAIL` | `--email` |
| `ncbi_api_key` | text | none | `PREPARE_NCBI_API_KEY`, then `NCBI_API_KEY` | `--api-key` |
| `temp_folder` | path (`~` expanded) | none | `PREPARE_TEMP_FOLDER` | `--temp-folder` |
| `log_file` | path (`~` expanded) | none | `PREPARE_LOG_FILE` | `--log-file` |
| `log_level` | level name | INFO | `PREPARE_LOG_LEVEL` | `--log-level`, `-q`, `-v` |
| `progress_every` | whole number, 0 or more | 50 | `PREPARE_PROGRESS_EVERY` | `--progress-every` |
| `log_host` | true or false | false | `PREPARE_LOG_HOST` | none |
| `min_free_gb` | number, 0 or more | 10 | `PREPARE_MIN_FREE_GB` | `--min-free-gb` |
| `prefetch_max_size` | whole number with optional K, M, G or T suffix | 100G | `PREPARE_PREFETCH_MAX_SIZE` | none |
| `assembly_memory` | `auto`, a fraction, or a size | auto | `PREPARE_ASSEMBLY_MEMORY` | `--assembly-memory` |
| `run_log` | true or false | true | `PREPARE_RUN_LOG` | none |

The flags exist on these commands:

- `--timeout`, `--temp-folder`: `download_sra` and `extract_target_reads`.
- `--lock-wait`: `download_sra` and `store_adopt`.
- `--email`, `--api-key`: `download_metadata`, `sra_info` and `validate_taxonomy`.
- `--min-free-gb`: `download_sra`. `--assembly-memory`: `extract_target_reads`.
- `--log-file`, `--log-level`, `-q`, `-v`, `--progress-every`: every command, before or after the command
  name.

A level name is DEBUG, INFO, WARNING, ERROR or CRITICAL, in any case. A true or false value may also be
written yes/no, on/off or 1/0.

What each setting does:

- `timeout`: how long `prefetch`, `fasterq-dump`, minimap2, samtools or megahit may run before it is
  stopped. A stopped download is classified as a network failure and retried. A `--version` probe
  always uses a fixed 30 s limit.
- `max_workers_cap`: the upper bound on the default number of parallel downloads, which is the CPUs
  available to the job divided by `--num-threads` (never more than 10, whatever the cap), and on the
  default number of datasets `sra_sketch` sketches at once. An explicit `--max-workers` is used as given.
- `lock_wait`: how long `download_sra` and `store_adopt` wait for another run's lock on the same
  accession before giving up on that accession.
- `registry_lock_wait`, `registry_lock_stale`: how long a command waits for the project registry lock,
  and how long a registry lock may go without a heartbeat before it is taken over. `registry_lock_stale`
  must be more than the registry heartbeat of 5 s.
- `dataset_lock_stale`, `lock_heartbeat`: the same two limits for the store dataset locks, the
  per-accession download locks of a plain project, and the index-build and extraction locks.
  `lock_heartbeat` must be less than `dataset_lock_stale`.
- `catalog_lock_wait`: how long a store catalogue write waits for the catalogue lock.
- `ncbi_email`, `ncbi_api_key`: sent with NCBI requests. A command that talks to NCBI stops with exit
  code 3 when no email address is set by any of the three routes. The key is never written to a log;
  it is shown as `(set)`.
- `temp_folder`: where `fasterq-dump` writes its temporary files, and where `extract_target_reads`
  writes its intermediate alignments. When it is not set, `download_sra` gives each accession a scratch
  folder of its own beside the download: `fastq/.prepare-tmp/<accession>_fqtmp` in a project without a
  store, `<data-root>/tmp/<accession>_fqtmp` with one. Set it to a local disk when the project or store
  is on a network filesystem and the node has local scratch space.
- `log_file`, `log_level`, `progress_every`, `log_host`: see [Logging](#logging) below.
  `progress_every` 0 turns the progress summaries off and logs one INFO line per item instead.
- `min_free_gb`: the free space a download of unknown size needs on each filesystem it writes to;
  0 turns the free-space check off. See [reads.md](reads.md).
- `prefetch_max_size`: the value passed to `prefetch --max-size`; a run whose `.sra` archive is larger is
  not fetched. A whole number of bytes, or a whole number with a K, M, G or T suffix (in either case).
- `assembly_memory`: the value passed to megahit `--memory`. `auto` gives 90% of the memory limit detected
  for the job (cgroup v2, cgroup v1, else `SLURM_MEM_PER_NODE` or `SLURM_MEM_PER_CPU`) and leaves
  megahit's own default when no limit is found. A size such as `32G` or `32000M` (binary units) or a whole
  number of bytes (at least 1M) is passed in bytes. A fraction such as `0.5` is passed as it is; megahit
  applies it to the whole node's memory, so under a scheduler give a size or `auto` instead.
- `run_log`: whether a command that keeps a run log adds a record of each run to
  `<project>/.prepare/runs/` (see [project_state.md](project_state.md)). Nothing is written in a folder
  without a project registry, and a run log that cannot be written is reported as a warning without
  changing the command's exit code. false turns the run log off.

## Order of precedence

A runtime setting is taken from the first of these that gives a value:

1. the flag, when the command has one and it was given;
2. the environment variable (`PREPARE_<NAME>`; for the API key, `PREPARE_NCBI_API_KEY` and then
   `NCBI_API_KEY`);
3. the `[runtime]` table of the config file;
4. the built-in default.

The root of the shared data store is taken from the first of these that names one:

1. `--data-root`;
2. the `PREPARE_DATA` environment variable;
3. `store.root` in the project registry (`prepare_registry.json`);
4. `data_root` in the `[store]` table of the config file.

A root found this way must hold the store marker (`prepare_store.json`), or the command stops with a
message naming the rule that produced the root. When none of the four names a root, the project works
without a store.

The two orders differ because the two kinds of value differ. A runtime setting tunes one run, and a
project has no record of its own for it, so the per-user config file is the last place before the
default. The store root describes where a project's data are, and a project that has joined a store
records that store in its registry; the registry therefore comes before the per-user config file, so
every command run in the project finds the store the project uses even when the user's default store is
a different one. The flag and `PREPARE_DATA` still come first, so a store that has been moved or is
mounted at another path can be named without editing the registry.

## Seeing the values in use

`prepare doctor` lists every runtime setting with its value and the source it came from (a flag, a
variable, `config [runtime] <key>`, or `default`), with the API key shown as `(set)`, and the store
root it resolved; the rule that produced the root is logged at INFO (`Resolved store root from
PREPARE_DATA: ...`). A command run at DEBUG (`-v` or `--log-level DEBUG`) logs the same list
at its start:

```
Runtime settings (value and source):
  subprocess_timeout = 0.0 (default)
  max_workers_cap = 4 (default)
  lock_wait = 3600.0 (--lock-wait)
  ...
  min_free_gb = 50.0 (config [runtime] min_free_gb)
```

The names in this list are the setting names of `SETTINGS`; they equal the config keys except for
`subprocess_timeout`, whose config key is `timeout`.

## Logging

Log lines go to stderr; a command's result (tables, JSON) goes to stdout. These options work on every
command, before or after the command name:

| Option | Effect |
|---|---|
| `--log-level LEVEL` | Console level: DEBUG, INFO (the default), WARNING, ERROR or CRITICAL |
| `-q`, `--quiet` | Only warnings and errors on the console (same as `--log-level WARNING`) |
| `-v`, `--verbose` | Debug lines and full tracebacks on the console (same as `--log-level DEBUG`) |
| `--log-file PATH` | Also append every line to `PATH`; its folder is created |
| `--progress-every N` | Log a progress summary every N items of a long run (default 50) |

`-q` and `-v` cannot be given together, and either one overrides `--log-level`. `store_status` lists every
dataset with `--details`; its former `--verbose` spelling still works but is deprecated. A flag
takes precedence over the environment variable and the `[runtime]` key of the same setting.

The log file is opened for appending, so several runs on one machine can share one file on a local disk.
On a network filesystem give each process its own file (for example one per SLURM array task, see
[hpc.md](hpc.md)), because appends from several hosts to one NFS file can interleave or be lost. The file
receives every line at INFO or above, also when the console is quiet, and DEBUG lines when the console is
at DEBUG. Each line names the host and process ID that wrote it:

```
2026-09-30 14:02:11 node17[48213] WARNING prepare.data.sra.retry: Failed to download SRR1234567: ...
```

When a command fails, the console shows the error on one line and the log file keeps the full traceback;
without a log file, rerun with `-v`. At DEBUG a run also logs the PREPARE version, its command line (with
the `--api-key` value hidden), the host, the process ID and, under SLURM, `SLURM_JOB_ID` and
`SLURM_ARRAY_TASK_ID`. `PREPARE_LOG_HOST=true` (or `log_host = true` in `[runtime]`) puts the host and
process ID on console lines too.

`download_sra`, `download_metadata` and `extract_target_reads` log a progress summary every
`--progress-every` items and at least every 5 minutes while items are finishing, then one closing line
with the totals and the time taken (for `download_sra`, after the retry pass):

```
download_sra: 150/2000 done (148 ok, 2 failed), 3.1/min, about 9 h 57 min left
download_sra: finished 2000/2000 (1990 ok, 10 failed) in 10 h 45 min
```

Lines about one accession or sample (a download starting, skipped or finished, its temporary folder, a
store link, a retry, an NCBI metadata request, an extracted sample) are logged at DEBUG;
`--progress-every 0` turns the summaries off and logs them at INFO instead. Warnings and errors about an
accession keep their own level. Some per-accession lines stay at INFO: the first "prefetch not found on
PATH" of a run and the lines about an unusual event (a download interrupted or redone with `--force`, a
truncated archive removed for the next attempt, a wait for running downloads to free space, a wait for a
lock another run holds).

## Checking the environment

`doctor` (group Environment in `prepare --help`) reports whether the environment can run PREPARE and
changes nothing:

```bash
prepare doctor                              # every check, one line each
prepare doctor --for extract_target_reads   # a tool this command needs is a failure, not a warning
prepare doctor --json                       # the same report as one JSON document
prepare doctor --network                    # also check that NCBI and Branchwater answer (10 s each)
```

The checks are:

- the Python and PREPARE versions;
- every external tool, with its path, version and the oldest supported version;
- that the config file parses, with each runtime setting and where its value came from;
- the shared data store, if one is configured, for its marker, write access and free space;
- free space at the project folder (`--project`, default the working directory), the temporary folder and
  the `.sra` cache against `min_free_gb`;
- the CPUs available against the node's total, the memory limit and any SLURM job variables;
- that the nearest project registry loads.

A missing tool is a warning unless `--for` names a command that needs it. A tool below its oldest
supported version, a config file or variable that does not parse, a store without its marker, an
unreadable registry and, with `--network`, an unreachable service are failures. `doctor` exits with 0 when
nothing failed (warnings allowed) and with 3 otherwise. It still runs when the config file is malformed,
which stops every other command, and reports the parse error as a failed check. `make doctor` runs it
from a checkout.

## Adding a setting

A new setting is one entry in `_SPECS` in `prepare/core/settings.py` (name, parser, default,
environment variable, description, and the flag's `dest` when a command has one) and one field in
`RuntimeSettings`. Code reads it with `settings.active().<name>`, or with `setting_for(args, name)`
inside a command so a flag in a hand-built namespace is still honoured. This page and the table above
list every entry of `SETTINGS`.
