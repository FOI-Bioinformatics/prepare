# Shared data store

This guide explains the user-level shared store: one copy of each downloaded SRA accession, linked into every project that uses it. Return to the [documentation index](README.md). Internal mechanics (locking, the catalogue, the journal, sidecar fields) are in [ARCHITECTURE.md](ARCHITECTURE.md).

A second organism project that studies the same metagenomes does not need its own copy of the reads. The store keeps one gzip-compressed copy of each accession, tracked in a small SQLite catalogue, and links it into each project.

## Finding the store

A store must exist before `download_sra` links a dataset or `store_adopt` runs. Without one, every command uses plain per-project folders. PREPARE looks for the store root in this order:

1. The `--data-root` flag.
2. The `PREPARE_DATA` environment variable.
3. `store.root` recorded in the project registry.
4. `[store] data_root` in `~/.config/prepare/config.toml` (`$XDG_CONFIG_HOME` is honoured).

If the configured root is unreachable, `status` and the analysis commands warn and continue without the store. `download_sra` and the `store_*` commands stop. `--data-root` is accepted by `download_sra`, `download_metadata`, `status`, `sra_profile`, `sra_report`, `sra_validate` and `extract_target_reads`. It never replaces `--fastq-folder`, which still names where the project expects its reads, as a folder or as the store's symlink.

## Commands

```bash
prepare store_init --data-root /data/prepare_store --set-default   # create a store, remember it
prepare store_adopt --fastq-folder fastq --dry-run                   # preview folding this project in
prepare store_adopt --fastq-folder fastq --move                      # move reads into the store, link back
prepare store_adopt --fastq-folder fastq --copy                      # copy into the store, keep the project folder
prepare store_status --json                                          # dataset and byte counts
prepare store_status --details                                       # also list every dataset in the store
prepare store_usage --accession SRR11011981                          # which projects used this run
prepare store_usage --project my_project                             # every dataset that project has used
prepare store_usage --organism GCF_000008025.1                       # every dataset used for that target genome
prepare store_usage --unused                                         # datasets in the store with no recorded usage
prepare store_gc --dry-run                                           # candidates for removal, nothing deleted
```

| Command | Purpose |
| --- | --- |
| `store_init` | Create a store at `--data-root` and record the project's use of it. `--set-default` also records the root as the user's default. |
| `store_adopt` | Fold an existing project `fastq/` folder into the store. |
| `store_status` | Dataset and byte counts. `--details` lists every dataset (the former `--verbose` is a deprecated alias). |
| `store_link`, `store_unlink` | Add or remove the project link of one accession. |
| `store_verify` | Check a dataset's files against the recorded size, md5 or spot count. |
| `store_reindex` | Rebuild the catalogue from the sidecar files on disk. |
| `store_usage` | Show which projects and genomes used an accession, or which datasets are unused. `--bytes-by-organism` totals bytes per genome. |
| `store_gc` | Report, or with `--yes` remove, datasets that nothing references. |

## Links in a project

Each linked accession appears as `fastq/<ACCESSION>`, a symlink to `<data-root>/sra/<ACCESSION>`. A dataset folder in the store may also hold `<ACCESSION>.sig.gz`, the sourmash signature `sra_sketch` cached for every project that links the dataset; `store_status` counts its bytes with the dataset. The link is relative when the store and the project share a parent folder and absolute otherwise. `--link-mode` (`auto`, `relative`, `absolute` or `copy`) overrides this choice. `store_init` and `store_adopt` add `fastq/` to the project's `.gitignore` when the project is a git repository.

With `--link-mode copy`, linking an accession that already has a copy in the project refreshes it. The new copy is made under a hidden name first, the earlier copy is moved aside only then, and it is put back if the swap does not happen, also on an interrupt. A relink therefore needs room for a second copy on the project's filesystem. The command checks that room first and fails with an error naming the accession when it is missing. A project folder that does not hold the store's `<ACCESSION>.json` is still refused.

With a store configured, `status` flags a wanted accession whose `fastq/<ACC>` links into the store but whose dataset is not yet `complete` or `unverified`. It reports such an accession as linked to a store dataset that is not complete, rather than as missing.

## Adopting existing downloads

`store_adopt --move` (the default) removes each accession's project folder and replaces it with a link to the store's copy. `store_adopt --copy` leaves the project's folder as it was, real and unlinked, whether the accession is new to the store or duplicates what the store holds. A `--copy` adoption still records the copying project as a user of the dataset, so `store_usage` and `store_gc` do not treat a copied dataset as unused.

A folder with no FASTQ files at all, such as an empty or interrupted download, is never adopted and is reported separately. `--metadata-folder` names the folder with `<ACCESSION>_metadata.xml` files, which supply the expected spot counts. `--lock-wait` bounds the wait for another project's lock. Adoption copies the data into the store before anything is removed from the project, so it briefly holds up to three copies of an accession and needs at least twice the folder's size free on the store's filesystem. The sidecar is written into the staging folder before it is published, so an adoption killed at any point leaves no store folder without a sidecar.

## Completeness verdicts

Every stored dataset carries a completeness verdict.

| Verdict | Meaning |
| --- | --- |
| `complete` | The read count per mate is at least 0.99 of NCBI's recorded spot count. |
| `partial` | Fewer reads than expected. `download_sra` does not use it unless `--accept-partial` is given, and downloads it again by default unless `--no-resume-partial` is given. |
| `unverified` | No expected spot count was found, because no metadata XML was present at download or adoption time. It is usable by default. |

The spot count comes from the lookup order described in [reads.md](reads.md), where the verdict rules are stated once. `store_verify --spots` and `status --reconcile` compute a verdict for a dataset that lacks one.

A store copy that is `partial` or `failed` is never linked into a project without `--accept-partial`. When `download_sra` meets such a copy:

- It refetches the copy unless `--no-resume-partial` is given, publishes the result to the store, and removes any project link into the store for that accession. A real project folder is never touched.
- It reports the accession as `incomplete: <ACC> store copy holds R of S spots; kept for --resume-partial; rerun with --accept-partial to use it`. This counts as a failed accession and the exit status is 1.
- A refetch that comes back with fewer reads than the published copy does not replace it, unless the refetch is verified `complete`.
- After two refetches in a row that gained no reads, which happens when NCBI's count cannot be reached, the copy is no longer fetched on every run. It is linked with `--accept-partial`, or refused with a message naming `--force`, which is then the only way to fetch it again. The counter is kept in the store sidecar and is shared by every project that uses the store.
- A result of `incomplete: ...`, or `partial in store; ...` under `--no-resume-partial`, is settled for the run. The retry pass does not fetch it again.

A copy recorded `unverified` whose sidecar read count is short against a spot count known since is treated as `partial` in the same way. `store_link` applies the same rule: it refuses these copies without `--accept-partial`, including an `unverified` copy whose read count is short against the project's spot count, and then records the sidecar as `partial`.

`status --reconcile` handles links into a store that is not mounted. When the `sra/` folder that the link points into does not exist, reconcile reports the link as store unavailable, logs one warning and leaves the link and the registry record unchanged. A dataset removed from a mounted store is still marked missing. A store-linked dataset that an earlier reconcile marked missing is recorded as linked again once its project entry is a symlink into the store.

## Verifying datasets

`store_verify` checks a dataset's files against its recorded size. `--md5` also compares each file's md5, and `--spots` compares the read count with NCBI's spot count. The spot count is read from the dataset's own sidecar, and otherwise from `<data-root>/metadata/<ACCESSION>_metadata.xml` and then from `metadata/<ACCESSION>_metadata.xml` in the calling project. A count found by `download_metadata`, either at download time or afterwards, is therefore picked up.

- **`--fix-state`** rewrites the sidecar and catalogue entry when a check finds a mismatch. It promotes a dataset to `complete` only after reading through its files, either through the `--spots` check or through a read-through that `--fix-state` performs itself. Matching size and md5 alone are not enough, since a file can still be a truncated or corrupt gzip stream.
- **`--rescan`** rebuilds a dataset's recorded file list from what is on disk before the check. Use it for files added or removed by hand.
- **`--json`** prints one document instead of the text table. It holds `root`, `checks` (which of `md5`, `spots`, `rescan` and `fix_state` ran), `counts` (the datasets per verdict, only verdicts that occurred), `datasets` (one row per accession checked, with its state before the check and, under `--fix-state`, a `fix` entry naming the action) and `fixed` (the accessions that `--fix-state` rewrote, with their state before and after).
- In text mode, `--fix-state` prints one line per rewritten dataset in the form `<accession>: fixed (<before> -> <after>)`.

## Rebuilding the catalogue

`store_reindex` rebuilds the SQLite catalogue from the sidecar files if it is lost. It replays the append-only journal under the store's `journal/` folder to restore which projects used which datasets.

- A store folder named like an SRA accession but without a sidecar does not stop the reindex. A sidecar is written from the file names and sizes alone, with state `failed` and an error naming `store_verify --rescan --fix-state <ACCESSION>`, and the dataset is catalogued.
- Such a dataset is not linked without `--accept-partial` until `store_verify --rescan --fix-state` has read its files and promoted it. A plain `--md5` check reports mismatches, because the rebuilt sidecar records no md5. When no project uses the dataset, `store_gc` offers it for removal.
- A folder whose lock is held, or whose name is not an SRA accession, is skipped with a warning. An unreadable sidecar still stops the reindex.
- If the replay restores no project while the rebuilt catalogue still holds datasets, `store_reindex` sets a `rebuilt_without_projects` flag, because every dataset would otherwise look unused. `store_gc` then refuses to run until each project that uses the store has run `store_init` or `store_link` again and `store_gc` is given `--accept-rebuilt`. A later `store_reindex` that restores at least one project on its own also clears the flag.

## Removing unused datasets

`store_gc` reports datasets that nothing references and, with `--yes`, removes them. Without `--yes` it only reports.

- It never removes a dataset that a project still links or that another run is working on.
- `--older-than DAYS` restricts it to datasets downloaded at least that many days ago.
- `--keep-partial` never removes a `partial` dataset.
- It also reports, and with `--yes` removes, leftover temporary artifacts under the store's `tmp/` folder from an interrupted download or adoption. This includes the bare `tmp/<ACCESSION>` staging folder of a killed download once that accession's lock is free.
- If the download's kept `.sra` archive or a temporary build folder grows past 1 GB combined, the download summary warns and names the folder, and `store_gc --dry-run` lists such leftovers as removal candidates.
- `--json` prints the report, or on a refusal `{"error": ...}`, as one JSON object. With no store configured, `store_gc`, `store_status` and `store_usage` print the same `{"error": ...}` shape under `--json`.
- On a store shared over a network filesystem or between machines, each project records the hostname it last ran on. A project not seen from the current machine looks stale even when it is active elsewhere. `store_gc` leaves a stale project's datasets alone unless `--include-stale` is given. The catalogue uses a rollback journal rather than WAL on such filesystems.

## Concurrent projects

A per-accession lock (`locks/<ACCESSION>.lock`, with a heartbeat) stops two projects from downloading the same accession into the store at once. A lock with no heartbeat for 10 minutes is treated as abandoned and taken over. `--lock-wait SECONDS` bounds how long `download_sra` and `store_adopt` wait for another project's lock. The default, 0, waits without a time limit for as long as the other project's heartbeat shows that it is still working. A positive value bounds the wait to that many seconds.

Without `--temp-folder`, `fasterq-dump` writes its scratch files to `<data-root>/tmp/<ACCESSION>_fqtmp` when a store is configured. See [hpc.md](hpc.md) for a store on a cluster.
