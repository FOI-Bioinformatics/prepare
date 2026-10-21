# Audit of the SRA download path, 2026-10-21

Scope: `prepare download_sra` from the command line to the files on disk, at commit 259df33
(v0.10.1): `cli/commands/sra.py`, `cli/commands/sra_verdicts.py`, `cli/commands/sra_record.py`,
the `data/sra/` package (`download.py`, `accession.py`, `retry.py`, `fastq.py`, `spots.py`,
`cleanup.py`, `space.py`, `run_report.py`, `store_handoff.py`), the store sidecar and
`utils/security.py`. The audit was read-only. Each finding was checked against the source by two
readers; none was reproduced on live data. Nothing in the repository was changed.

The second half of this document describes the path as built, because no other document does so
at this level.

## 1. Findings

Severity: H = data loss, wrong result recorded as correct, or a wrong exit code; M = wrong
behaviour without an error, or a misleading error; L = cosmetic, dead code, or duplication.

### High

| # | Where | Defect | Failure scenario |
|---|-------|--------|------------------|
| H1 | `data/sra/accession.py:242-284`; `data/sra/store_handoff.py:548-550` | A truncated download (mate-1 read count below 99% of `run_total_spots`) returns `True` from `_handle_download_output` in a plain project. The folder is published, the registry records `downloaded`, the accession is not retried and the exit code is 0. The same data through the store path ends as `(False, "incomplete: ...")`, exit code 1. | A run with `run_total_spots` 1,000,000 whose fasterq-dump stops at 500,000 reads: in a plain project the run reports success and the half dataset enters every later analysis; with a store it is a failure. |
| H2 | `data/sra/retry.py:133-137, 455-462` | A worker that raises is stored as `(acc, None)`. `_process_download_results` then fails to unpack it and records the message `Error: cannot unpack non-iterable NoneType object`. That text is what the report CSV, `download_run.json`, `failure_reason`, the exit-code choice and the retry split see; only the registry callback gets the real message. `tests/test_data_sra_retry.py:69-72` pins the placeholder text. | A `NetworkError` from one worker is classified `unknown` instead of `network`, so the exit code is 1 instead of 4 and the retry pass treats it as a generic failure. |
| H3 | `cli/commands/sra_verdicts.py:105-116`; `data/sra/store_handoff.py:84-94` | `--no-verify-downloads` removes the spot counts from the registry inputs, but `_resolve_expected_spots` falls back to the store sidecar and the metadata XML, so a store download is always verified. | A user who passes `--no-verify-downloads` because NCBI's spot count is known to be wrong for a run still gets `incomplete: ...` and no link in the project. |
| H4 | `data/sra/accession.py:242` | The plain path computes the verdict on the uncompressed files and never reads the `.gz` files again or records a digest. The store path reads each compressed file once in `build_sidecar`. | A pigz failure that leaves a short `.gz` is logged as a warning and the folder is published as complete. |

### Medium

| # | Where | Defect | Failure scenario |
|---|-------|--------|------------------|
| M1 | `data/sra/accession.py:248-262, 281-282`; `store_handoff.py:515` | An interrupt during compression returns success with a mixed folder of `.fastq` and `.fastq.gz` files. The store sidecar then labels compression `gzip` when any one file is gzipped. | Ctrl-C while pigz runs on mate 2 publishes `SRR1_1.fastq.gz` beside `SRR1_2.fastq`; a later step that expects one layout per dataset reads half of it. |
| M2 | `data/sra/accession.py:310-322` | `_cached_sra_archive` globs `<ACC>.sra*`. A prefetch leftover (`.sra.lock`, `.sra.prf`, `.sra.tmp`) sorts before `.sralite` and is handed to fasterq-dump as its source. | The path does not end in `.sra`, so `run_secure` rejects it as an accession and the run fails with `Invalid SRA accession format` for a valid accession. |
| M3 | `data/sra/accession.py:449-455, 321` | When prefetch exits 0 without writing an archive, fasterq-dump is run on a path that does not exist and the failure is classified from its message. | A run above `--max-size` is reported as a fasterq-dump failure, not as a size refusal. |
| M4 | `data/sra/retry.py:245-283, 298` | The retry pass runs on the main thread and has no interrupt collection. Queued workers' `"interrupted"` results are recorded as `failed` with an attempt counted, and the back-off `time.sleep(2**retry)` is not interruptible. | Ctrl-C during the second retry round records attempts that never ran and waits out the sleep. |
| M5 | `data/sra/accession.py:488-503, 219-220, 566` | A failed attempt keeps `.prepare-tmp/<ACC>_temp` for inspection and nothing removes it until the next attempt on the same accession. | Ten failed accessions leave ten scratch folders of partial FASTQ; the only signal is the warning when transient folders exceed 1 GiB. |
| M6 | `data/registry.py:357-361` | A failed forced redownload keeps the old folder on disk but records `files=[]` and `bytes_total=0`. | `status` reports the dataset as failed with no files while the previous complete copy is still there. |
| M7 | `data/sra/store_handoff.py:552, 636` | When the project folder `fastq/<ACC>` is a real folder, `link_dataset` refuses it only after the store copy is published and catalogued. The error is classified `unknown` and retried. | A project that downloaded before a store existed re-fetches the run into the store, then fails twice on the link. |
| M8 | `cli/commands/sra.py:423-435`; `utils/security.py:330` | An explicit `--max-workers` is not bounded by `MAX_CONCURRENT_DOWNLOADS` or `max_workers_cap`; `--threads` only has to be digits, so `0` passes. | `--max-workers 64` starts 64 fasterq-dump processes; `--num-threads 0` reaches fasterq-dump. |

### Low

| # | Where | Defect |
|---|-------|--------|
| L1 | `data/sra/retry.py:311-312` | `failed_accessions.txt` from an earlier run is not removed by a run with no failures. |
| L2 | `data/sra/retry.py:319-323`; `cli/commands/sra.py:287-291` | The hint "To retry only failed accessions" is logged twice. |
| L3 | `cli/commands/sra.py:152-158` | `--verify-downloads` is `store_true` with `default=True`, so the flag does nothing. |
| L4 | `data/sra/accession.py:421-431`; `data/sra/download.py:383` | `redownload_truncated` is applied to every store accession, so the archive cache is wiped for accessions that are not truncated. |
| L5 | `core/constants.py:33`; `utils/security.py:57`; `store_handoff.py:31`; `retry.py:252, 378, 565` | Unused: `DEFAULT_MAX_RETRIES`, `FASTERQ_DUMP_BOOLEAN_FLAGS`, `MAX_SUBPROCESS_TIMEOUT`, `_publish_store_dataset`, and the `downloader or download_accession` fallback, which only tests reach. |
| L6 | `store_handoff.py:194` vs `sra_verdicts.py:258`; `store_handoff.py` `_metadata_xml`/`_xml_folders` vs `spots.spots_from_xml`; `sra.py:522` vs `store_handoff.py:516`; `fastq.count_fastq_reads` vs `_NewlineCounter` | The same logic exists twice with different type rules, and the store metadata folder is searched twice. |
| L7 | `store_handoff.py:523` | `fasterq_dump_version()` runs `fasterq-dump --version` for every store download. |
| L8 | `utils/security.py:421-439`; `accession.py:300, 303` | `--progress` output of prefetch and fasterq-dump is captured into memory; the flag is given in direct mode only. |
| L9 | `data/sra/accession.py:28` | `_PREFETCH_FALLBACK_LOGGED` is process-wide, so "the run's first accession" means the process's first. |

### Test gaps

The fake tools in `tests/helpers_processes.py:101-128` never fail, never write partial output and
never write an unpaired file. The following paths have no end-to-end test:

- a tool failure through classification, retry and the registry record;
- a truncated verdict, in either a plain project or a store;
- `--timeout` to `SecurityError` to the network classification and a retry;
- prefetch leftovers beside the archive;
- an interrupt inside the retry pass;
- `--no-verify-downloads` with a store;
- compression interrupted, then published, then recorded.

## 2. The path as built

### 2.1 Steps for one accession

From the command (`cli/commands/sra.py:627`): `require_tools(["fasterq-dump"])`, the worker
count, the registry and store, then `registry_inputs` (`sra_verdicts.py:83-123`) returns the
excluded accessions, the expected spot counts, the truncated set and the run sizes.
`_download_batched` (`sra.py:526`) opens a `registry_batch` and the termination context and calls
`download_sra` (`data/sra/download.py:219`), which reads the accessions, the blacklist and the
existing downloads, applies `--max-downloads`, picks the worker (`store_handoff._store_downloader`
with a store, `accession._project_download` without) and runs `retry._download_with_retries`.

Plain project (`accession.py:575`):

1. Take `<fastq>/.locks/<ACC>.lock`; check the folder again under the lock.
2. Wipe `<fastq>/.prepare-tmp/<ACC>_fqtmp` (the fasterq-dump scratch folder) unless
   `--temp-folder` is given.
3. `_stage_and_publish` wipes `.prepare-tmp/<ACC>` and calls `download_accession` with that
   staging folder.
4. `download_accession` (`accession.py:365`) removes an existing folder when forcing, removes the
   archive cache on a redownload, recreates `<staging>/<ACC>_temp`, runs prefetch when requested
   and on PATH (`prefetch -O <cache> --max-size <prefetch_max_size> --progress <ACC>`), runs
   fasterq-dump (`--split-3 --skip-technical --threads N -O <ACC>_temp [--temp T] <archive>` with
   prefetch, or `--threads N --progress <ACC> -O <ACC>_temp [--temp T] --split-3 --skip-technical`
   without), then `_handle_download_output` moves the files into `<staging>/<ACC>/`, computes the
   verdict on the uncompressed files and compresses each with pigz (`-p N -f`) or Python gzip.
   On success the archive cache is removed unless `--keep-sra`.
5. `verify_held(lock)`, then `publish_folder` renames the old folder to `<tmp>/<ACC>_old`, the
   staged folder into place, and deletes the old one.

Store (`store_handoff.py:557`): resolve the expected spot count; without `--force` run the
precheck lock-free and again under `dataset_lock` (a ready copy is linked; an incomplete copy is
refused or linked according to `--resume-partial` and `--accept-partial`; two refetches that
gained nothing settle the dataset). `_store_fetch` forces a full refetch whenever `sra/<ACC>`
exists, with cache `store.tmp/.sra-cache` and scratch `store.tmp/<ACC>_fqtmp`, and calls
`download_accession` into `store.tmp`. `build_sidecar` reads each file once for bytes, md5 and
records and sets the state; after `verify_held` and `_publish_decision` the copy is published to
`sra/<ACC>`, catalogued, and linked into the project unless incomplete.

### 2.2 Layout and names

- Plain project: `<fastq>/<ACC>/<ACC>_1.fastq.gz`, `<ACC>_2.fastq.gz`, and the unpaired
  `<ACC>.fastq.gz` that `--split-3` writes.
- Store: `<store>/sra/<ACC>/` with the same files plus `<ACC>.json`; the project's `fastq/<ACC>`
  is a symlink or a copy.
- Transient names (`cleanup.py:21-31`): `<ACC>_temp`, `<ACC>_fqtmp`, `.sra-cache`,
  `.prepare-tmp`, `<tmp>/<ACC>_old`.
- Run files: `<fastq>/failed_accessions.txt` and `<fastq>/download_run.json`.

### 2.3 Completeness

`verify_download` (`fastq.py:325`) counts newlines divided by four in the mate-1 file (or the
single-end file) plus the unpaired file when it sits beside a pair. `verdict_for_count`
(`spots.py:44`) gives `complete` at a ratio of 0.99 or more, `truncated` below, `unverified` when
no spot count is known. The spot count is looked up in the registry `run_total_spots`, then the
sidecar, then `<ACC>_metadata.xml`, then the verdict already on file. md5 and size are computed
only in the store path and never compared with NCBI's `Run_MD5`, which is only stored.

### 2.4 Resume, retry, interrupt

No partial state is resumed: the staged folder, `_temp` and the scratch folder are wiped before
every attempt. The prefetch archive is reused on a later non-forced attempt and prefetch resumes
it itself. `--max-retries` defaults to 1; the retry pass is sequential on the main thread with a
back-off of `2**retry` seconds. Not-found, insufficient-space and settled store outcomes are
not retried. A disk-full error aborts the run. On interrupt the first pass sets `stop`, cancels
pending futures, terminates the children (SIGTERM, SIGKILL after 5 s), collects finished results
except `"interrupted"` ones and re-raises; the CLI exits 130 and still writes the summary.

| Setting | Default | Controls |
|---|---|---|
| `subprocess_timeout` (`--timeout`) | 0 | prefetch, fasterq-dump and pigz; a timeout is classified as network |
| `lock_wait` (`--lock-wait`) | 0 | wait for a held lock while the holder is alive |
| `dataset_lock_stale` | 600 s | when a store lock without heartbeat is taken over |
| `lock_heartbeat` | 10 s | heartbeat interval |
| `max_workers_cap` | 4 | upper bound of the default worker count |

### 2.5 Parallelism

One `ThreadPoolExecutor` runs the first pass; each worker runs its own fasterq-dump with
`--threads N` and pigz with `-p N`. The default worker count is
`min(MAX_CONCURRENT_DOWNLOADS, cpus // num_threads, max_workers_cap)`. The stop token travels
from `BaseCommand.run` (`args._termination`) through `download_sra`, the worker, the lock
helpers, `download_accession`, `run_secure(stop=)` and `compress_fastq(stop=)`.

### 2.6 Records

Per result, `record_download` writes `attempts`, `state`, `date`, `message`, `bytes_total`,
`files[{path, bytes, mtime}]` (only when downloaded), `complete{verdict, ratio, expected_spots,
reads_r1}`, `source` and `store_name`; `set_download_timing` adds `started` and `seconds`. After
the run, present, blacklisted and skipped accessions are recorded. The run log summary holds
`accessions`, `downloaded`, `failed`, `already_present`, `blacklisted`, `skipped`, `attempts`,
`failures_by_reason` and `aborted`; the detail holds `failed{acc:{reason, attempts, message}}`
and `downloaded{acc:{status, seconds, attempts}}`.

### 2.7 Allowlist

| Tool | Allowed flags | Value checks |
|---|---|---|
| `fasterq-dump` | `--threads --progress -O --temp --split-files --split-3 --skip-technical --include-technical --force --version` | `--threads` digits only; `-O`, `--temp` path-checked |
| `prefetch` | `-O --max-size --progress --resume --version` | `--max-size` shell characters only |
| `pigz` | `-p -f -k --version` | file argument not path-checked |

Positional arguments must be a `.sra`/`.sralite` path or match `^[A-Z]{3}[0-9]+$`. `fastq-dump`,
`vdb-dump`, `vdb-validate` and `sam-dump` are not allowed and not referenced. `TOOLS` lists
`fasterq-dump` (required), `prefetch` and `pigz` (optional), all sra-tools 3.0 or later.

## 3. Fix plan

The fixes are planned in `docs/superpowers/plans/2026-10-21-sra-download-fixes.md`, in four
waves: verdict consistency (H1, H3, H4, M1), error reporting (H2, M4, M7), robustness (M2, M3,
M5, M6, M8, L1-L4) and dead code with the missing process tests (L5-L7).
