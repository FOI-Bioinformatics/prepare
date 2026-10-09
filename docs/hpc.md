# Running PREPARE on a cluster

This page describes how to run the long steps of PREPARE, `download_sra` and `extract_target_reads`,
as SLURM batch jobs: checking a node before the first job, splitting a download over an array job,
what happens at the walltime limit, how exit codes drive a resubmission, and what to watch for when
the project or the shared data store lives on a network filesystem. PREPARE does not submit jobs
itself; the scripts below are templates to adapt to the local cluster.

The settings named here (`--timeout`, `--lock-wait`, `--min-free-gb`, `--assembly-memory`, and the
`PREPARE_*` variables) are described in [configuration.md](configuration.md).

## Before the first job

Run `doctor` inside an allocation of the same size as the jobs, so it sees the job's CPUs, memory
limit and temporary folder rather than those of the login node:

```bash
srun --cpus-per-task=8 --mem=16G --time=00:10:00 prepare doctor --for download_sra
srun --cpus-per-task=16 --mem=64G --time=00:10:00 prepare doctor --for extract_target_reads
```

`--for` turns a tool the command needs into a failure rather than a warning, and `doctor` then exits
with 3. It also reports the config file and every runtime setting with its source, the shared data
store (marker, write access, free space), free space at the project, the temporary folder and the
`.sra` cache, the CPUs available against the node's total, the memory limit, and the SLURM variables
it finds. `--json` writes the same report as one JSON document.

## Downloads as an array job

Split the accession list into one shard per array task and give each task its own log file and
report file. `split -n l/N` makes N shards without breaking a line:

```bash
N=10
mkdir -p shards logs reports
split -n l/$N -d accessions.txt shards/shard_      # shards/shard_00 ... shards/shard_09
sbatch download.sbatch
```

`download.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=mq-download
#SBATCH --array=0-9
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --signal=B:TERM@300
#SBATCH --output=logs/slurm_%A_%a.out

shard=$(printf 'shards/shard_%02d' "$SLURM_ARRAY_TASK_ID")

exec prepare --log-file logs/dl_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log download_sra \
    --accessions-file "$shard" --num-threads 4 --max-workers 2 \
    --temp-folder "$TMPDIR" --timeout 0 --lock-wait 3600 \
    --report-file reports/dl_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.csv
```

Notes on the template:

- The array range is `0` to `N-1`. Append `%K` (`--array=0-9%4`) to run at most K tasks at once, which
  limits the load on NCBI and on the filesystem.
- `--num-threads 4 --max-workers 2` uses the 8 CPUs of the allocation: two accessions at a time, each
  `fasterq-dump` with four threads. Without `--max-workers` the number is the CPUs available to the job
  divided by `--num-threads`, at most 4 (`PREPARE_MAX_WORKERS_CAP`).
- `--temp-folder "$TMPDIR"` puts the `fasterq-dump` temporary files on the node's local disk, when the
  site sets `TMPDIR` to one. Check that it is set, and large enough: the free-space check asks for about
  8 times an accession's `.sra` size there while it converts. Without `--temp-folder` (or
  `PREPARE_TEMP_FOLDER`) the scratch files go to `fastq/.prepare-tmp/<ACCESSION>_fqtmp` in a project
  without a store, and to `<data-root>/tmp/<ACCESSION>_fqtmp` with one, that is, onto the project's or the
  store's filesystem; when that is a network filesystem, point the scratch at local disk as here. The
  free-space check and `doctor` measure whichever of these folders `fasterq-dump` writes to.
- `--timeout 0` lets a tool run for as long as it needs; the walltime is the limit. This is the default
  and is written out here only to make it visible.
- `--lock-wait 3600` gives up on an accession after one hour when another run holds its lock (a second
  project downloading the same run into the shared store, for example). The accession is recorded as
  failed with a `locked:` message and is picked up by the next run.
- SLURM does not create the folder of `--output`, so `logs/` must exist before `sbatch`. `--log-file`
  creates its own folder, `--report-file` does not.
- Every task writes to the same project registry, `prepare_registry.json`, under its lock. A task
  queues its outcomes and writes them in batches, so the lock is held briefly. On a busy network
  filesystem, raise the wait for that lock with `export PREPARE_REGISTRY_LOCK_WAIT=120` in the script
  (default 30 s).

Accessions already present are skipped, so a shard can be run again at any time: a rerun downloads only
what is missing.

## Signals and the walltime limit

`--signal=B:TERM@300` asks SLURM to send `SIGTERM` to the batch shell 300 seconds before the time limit.
`exec` replaces the shell with `prepare`, so the signal reaches PREPARE itself; without `exec` the
shell receives it and PREPARE runs on until SLURM kills the job.

On the first `SIGTERM` (or `SIGINT`, `SIGHUP`) `download_sra`:

- starts no further accession, and stops the `prefetch` and `fasterq-dump` processes that are running;
- writes the outcome of every accession that finished to the project registry;
- does not publish the accession that was cut off: its files are built in a staging folder and moved
  into place only when complete, so it is not reported as present and the next run downloads it again;
- exits with status 130.

It still writes `--report-file` and `fastq/download_run.json` for an interrupted run, listing the
accessions that reported a result before the signal (`aborted` is `"interrupted"`, exit code 130). It
does not write `fastq/failed_accessions.txt` for an interrupted run; the registry and the log file hold
what was done. A second signal is logged and does not interrupt the registry write; a third stops at once.

`extract_target_reads` records each sample in the registry as it finishes. On a signal it stops the
running tool, keeps every sample recorded so far and exits with 130. The sample that was cut off is done
again by the next run.

Each tool runs in a process group of its own, which PREPARE signals as a whole on a timeout or a
termination, so a process the tool started (megahit's `megahit_core`) stops with it. A `SIGKILL` sent
only to PREPARE's process group does not reach a running tool. SLURM ends every process of a job
through its process tracking (`proctrack/cgroup` on most clusters); with `proctrack/pgid` a tool can
outlive a job that was killed rather than signalled.

To resume after the time limit, submit the same script again. If a log file does not end with
`Download interrupted by the user`, the job was killed before the final registry write finished;
give it more time with a larger value than 300.

## Exit codes and job dependencies

| Code | Meaning | In a pipeline |
|---|---|---|
| 0 | Success | continue |
| 1 | Failure: bad input, a dataset not found, a tool error | stop and read the log |
| 2 | Usage error or renamed command | fix the script |
| 3 | Configuration: a missing or too old tool, a malformed config file, no NCBI email | fix the environment |
| 4 | Retryable: a network failure, or a wait for a registry, catalogue or store lock that gave up | resubmit later |
| 130 | Interrupted by a signal, including the walltime signal | resubmit to resume |

`download_sra` exits with 4 only when every accession that failed did so for a network reason (a
connection or timeout error, or a tool stopped by `--timeout`). If any accession failed for another
reason (not found, disk full, not enough free space, locked by another run, including one given up
after `--lock-wait`) it exits with 1.

Chain the steps with `--dependency=afterok`, which starts the next job only when every task of the
previous one exited with 0:

```bash
dl=$(sbatch --parsable download.sbatch)
sbatch --dependency=afterok:$dl extract.sbatch
```

To resubmit the download while it ends with 4, run a loop on the login node (or in a small job of its
own). `sbatch --wait` returns when the job ends, with the job's exit code; for an array job this is the
highest exit code of its tasks.

```bash
#!/bin/bash
# resubmit_download.sh: rerun the download array while it ends with exit code 4
for attempt in 1 2 3 4 5; do
    sbatch --wait download.sbatch
    rc=$?
    if [ "$rc" -ne 4 ]; then
        break
    fi
    echo "download attempt $attempt ended with 4; resubmitting in 30 minutes" >&2
    sleep 1800
done
exit "$rc"
```

Because the highest code wins, a task that ends with 130 (walltime) or with 1 alongside one that ends
with 4 gives 130 or 4 respectively; read the per-task logs before deciding. Each resubmission downloads
only the accessions still missing, so a rerun costs little.

## A shared data store on a group filesystem

A shared store (see "Shared data store" in the README) keeps one copy of each downloaded run for every
project of a group. On a cluster it normally lives on a group filesystem:

```bash
prepare store_init --data-root /proj/group/prepare-store --set-default
export PREPARE_DATA=/proj/group/prepare-store     # in the job script, or rely on the registry
```

The store root is found from `--data-root`, then `PREPARE_DATA`, then `store.root` in the project
registry (written by `store_init`), then `[store] data_root` in the user config file, so a project that
ran `store_init` finds its store in every job without further settings.

- Every member of the group needs write access to the whole store: a waiting run may take over a lock
  file that another user's run left behind. Make the root group-writable with the set-group-ID bit
  (`chmod -R g+rwX /proj/group/prepare-store` and `chmod g+s` on its folders) and set `umask 0002` in
  the job scripts.
- Each dataset has a lock, `<data-root>/locks/<ACCESSION>.lock`. Two jobs, of one project or of two,
  that want the same accession do not both download it: the second waits (up to `--lock-wait`) and then
  finds it present. A lock whose heartbeat has stopped for 10 minutes (`PREPARE_DATASET_LOCK_STALE`) is
  taken over. A holder that died on the same host is taken over at once; one on another host cannot be
  checked, so its lock is taken over only after the full 10 minutes.
- With a store, the finished FASTQ files are staged under `<data-root>/tmp/` and moved into place with
  one rename, which needs the staging folder on the same filesystem as the store. Keep `--temp-folder` on
  local disk; it holds only the `fasterq-dump` scratch files.
- `store_gc` keeps a dataset for one day after it was last linked or downloaded, and never removes one a
  registered project still uses.

## Network filesystems

PREPARE's locks are plain files created with `O_EXCL`, and the store catalogue is an SQLite database
in `DELETE` journal mode behind its own lock file. This works on network filesystems, with these
caveats:

- **O_EXCL**: file creation with `O_EXCL` is atomic on NFSv4 and on NFSv3 with a current Linux client.
  It is not atomic on NFSv2; do not put a project or a store on such a mount.
- **Clock skew**: a waiter judges a lock stale by comparing its own clock with the lock file's
  modification time, which the file server may set from its own clock. Keep the nodes and the file server
  synchronised (NTP). A skew larger than a stale window (120 s for the registry, 600 s for datasets) can
  make a live lock look stale, or delay the takeover of a dead one by the size of the skew. NFS clients
  also cache file attributes for up to about a minute by default; the stale windows are chosen above
  that, and should not be lowered much below the defaults on NFS.
- **Temporary files**: do not point `--temp-folder` (or `PREPARE_TEMP_FOLDER`) at a network
  filesystem. `fasterq-dump` writes and rereads its scratch files many times; on NFS this is slow and
  loads the file server. Use the node's local disk (`$TMPDIR`).
- **Log files**: give each process its own `--log-file`, as the array template does. Appends from
  several hosts to one file on NFS are not atomic, and lines can interleave or be lost. Several
  processes on one host may share a log file on a local disk.
- **SQLite**: the catalogue does not use WAL mode, which SQLite documents as unsafe on network
  filesystems; writes are serialised by `catalog.sqlite.lock` and wait up to 60 s
  (`PREPARE_CATALOG_LOCK_WAIT`).

## Targeted extraction and assembly

Run one `extract_target_reads` job per target genome, for example as an array over a list of genome
IDs. Within one job the samples are processed one after another; the index for the genome is built once
and reused. Two jobs on the same genome do not map a sample twice (a sample another job is working on is
reported as skipped), but the skipped samples are then not done by that job and need a later rerun, so
splitting one genome over several jobs gains little.

`extract.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=mq-extract
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --signal=B:TERM@300
#SBATCH --output=logs/slurm_%A_%a.out

genome=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" genomes.txt)

exec prepare --log-file logs/ex_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log extract_target_reads \
    --parsed-containment parsed_containment.txt --genome-id "$genome" --genome-fasta "genomes/$genome.fna" \
    --threads "$SLURM_CPUS_PER_TASK" --temp-folder "$TMPDIR" --assemble --assembly-memory auto
```

`--assembly-memory auto` (the default) passes megahit `--memory` as 90% of the memory limit of the job,
read from the job's cgroup (v2, then v1) or from `SLURM_MEM_PER_NODE` (`--mem`) or
`SLURM_MEM_PER_CPU` (`--mem-per-cpu`). Without it megahit sizes itself from the
node's total memory and can be killed for exceeding the job's `--mem`. A fixed size (`--assembly-memory
56G`) also works; a fraction (`0.9`) does not help under a scheduler, since megahit applies it to the
whole node. megahit uses `--threads` threads on Linux unless `--assembly-threads` is given.

A sample already extracted or assembled with the same genome, preset and threshold is skipped, so the
job can be resubmitted after the walltime limit to finish the remaining samples.

Assemblies resume safely after a walltime signal or a kill. megahit writes into a hidden staging folder
beside `<genome_id>_assembly`, and the finished result replaces the previous folder by rename, so an
interrupted run never leaves a half-written assembly in place of a complete one. A staging folder left by a
killed job is removed when the next job reaches that sample, and a folder without `final.contigs.fa` is
assembled again. Each finished assembly carries a marker, `.prepare-assembly.json`, recording the reads
and settings it was built from; a resubmitted job reuses an assembly whose marker matches and redoes one
whose reads or settings have changed (see "Targeted Read Extraction Before Assembly" in the README). A
sample whose assembly fails is reported and the job continues with the next one, then exits with 1, so an
`afterok` dependency on the job does not start; a sample another job is working on is skipped and does not
fail the job. Each sample's assembly step loads the project registry once, also when its assembly is
reused, which costs about one second per sample on a registry of about 15,000 datasets; a resubmitted job
in which every assembly is reused spends that time and little else.

## Reconciling after an upgrade to 0.8.0

From 0.8.0, `status --reconcile` re-checks every download verdict recorded as `unverified` once a spot
count is known for it. A store link costs nothing for this, since the read count comes from its sidecar,
but a plain project download is counted again from its mate-1 file and its file of unpaired reads. The
first reconcile after upgrading therefore reads, once, the mate-1 file of every `unverified` plain download
that has a spot count; later runs read nothing again. On a large project, run that first reconcile as a job
rather than on a login node. A project link into a store that is not mounted on the node is reported as
store unavailable and left unchanged, so a reconcile on a node without the store mount does not mark its
datasets missing. The same first reconcile drops every assembly record dated before its extraction record,
which a 0.7.0 registry holds for a sample extracted again after its assembly, and lists it under
`assemblies_dropped`; a later `extract_target_reads --assemble` job builds those assemblies again.
