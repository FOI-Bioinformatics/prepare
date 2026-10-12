# Audit fixes 0.10.1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix every confirmed and likely bug in `docs/superpowers/audits/2026-10-12-code-audit.md`, close the CI and packaging gaps, and do the structural work the audit recommends, released as 0.10.1.

**Architecture:** Six waves in the audit's fix order. Wave 1 closes security and data-loss holes by validating accessions at every boundary and deriving project folders from the registry location instead of the cwd. Wave 2 corrects wrong results with value-checked tests and replaces the three tests that pin bugs. Wave 3 carries `force` and the retry prefix through the download result path. Waves 4 and 5 fix reports and CLI behaviour and bring CI to parity with `make check`. Wave 6 is structural: one NCBI pacer, a finished plugin base, module splits, and the CSS prefix rename. Each task is one commit on a branch `fix/audit-0.10.1`; the branch is squash-reviewed at the end.

**Tech Stack:** Python 3.12, pandas, numpy, scipy, pytest. Run everything with `conda run -n metaquest <cmd>` from `/Users/andreassjodin/Code/prepare` (the shell cwd resets to the old folder, so prefix every command with `cd /Users/andreassjodin/Code/prepare &&`).

**Spec:** `docs/superpowers/audits/2026-10-12-code-audit.md`. Bug IDs (H1-H5, M1-M16, L1-L11, "likely") refer to that file.

## Context

The 2026-10-12 audit of v0.10.0 found 32 confirmed and about 25 likely bugs. Five are high severity: a path traversal that deletes folders outside the project from a bad accessions line, a `store_gc --include-stale` data-loss path caused by an unresolved registry path, a wrong PERMANOVA statistic, registry corruption from `status --reconcile` run in a subfolder, and metadata colouring that never worked in PCA/t-SNE. The user chose to fix all of it (all six audit steps, every likely bug) under version 0.10.1.

## Global Constraints

- ASCII only in `prepare/`, `tests/`, `scripts/`, `Makefile`, `setup.cfg`, `pyproject.toml` (`scripts/check_ascii.sh`).
- No `print(` outside `prepare/cli/base.py`; commands use `self.emit`.
- Every output file is written through `write_text_atomic` / `atomic_path` / `open_atomic` (`scripts/check_atomic_writes.sh`, allowlist `scripts/atomic_writes_allowlist.txt`).
- Modules stay under 800 lines; cyclomatic rank below D (`make check`).
- black at 120 columns, flake8, mypy clean.
- Commit trailer: `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`. Never `git add -A` (AGENTS.md is untracked on purpose).
- Modest scientific language in docs and docstrings.
- Public CLI flags and output file names do not change in this release. Behaviour changes are listed in CHANGELOG under 0.10.1.
- Gate before every commit: `conda run -n metaquest make check` and the owned tests; full `make test` at the end of each wave.

## Review Focus

1. An accession string containing `/`, `..`, or an absolute path, at every entry point: must be refused with `ValidationError` before any filesystem action (Tasks 1.1, 1.2).
2. Any command run from a subfolder of the project: must act on the project found by walking up, never on `<cwd>/fastq` (Task 1.4).
3. A retry result (`"Retry 1: linked from store"`) must be recorded exactly like the unprefixed result (Task 3.2).
4. A metadata table whose index is not in frame order, or has duplicate rows: PCA colouring must align by sample name and the explorer must not crash (Tasks 4.1, 4.3).
5. A statistic on degenerate input (one group, a NaN label, zero-length read, empty list): must return a defined value or a clear error, never a traceback (Tasks 2.1, 3.4, 2.6).

## File Structure

New files:
- `prepare/data/accessions.py`: the one accessions-file reader and `require_accession()`.
- `prepare/data/ncbi_pacing.py`: `RequestPacer` shared by the three NCBI clients (wave 6).
- `prepare/data/registry_io.py`, `prepare/data/registry_scan.py`: split out of `registry.py` (wave 6).
- `prepare/data/read_extraction_export.py`: FASTQ export and coverage, split out of `read_extraction.py` (wave 6).
- `prepare/cli/commands/sra_record.py`: the result recorder, split out of `cli/commands/sra.py` (wave 6).
- `prepare/py.typed`: empty marker.
- `scripts/check_version.py`: tag/version consistency.
- `tests/test_accessions.py`, `tests/test_audit_regressions.py` (one class per wave for cross-cutting regressions), plus additions to existing test modules named in each task.

---

## Wave 1: security and data loss

### Task 1.1: One validated accessions reader (H1, likely: comment lines, duplicates)

**Files:**
- Create: `prepare/data/accessions.py`
- Modify: `prepare/cli/base.py:26-48`, `prepare/data/sra/download.py:322-326`, `prepare/data/metadata.py:288` (`_read_accessions_file`), `prepare/data/sra/accession.py:355-405` (top of `download_accession`)
- Test: `tests/test_accessions.py`, `tests/test_data_sra.py`

**Interfaces:**
- Produces: `prepare.data.accessions.read_accessions_file(path, validate=True) -> List[str]` (strips, skips blank and `#` lines, de-duplicates in order, raises `ValidationError` naming the file, line number and offending text when `validate` and a line does not `re.fullmatch(SRA_ACCESSION_PATTERN)`; raises `DataAccessError` on OSError). `require_accession(accession: str) -> str` (fullmatch or `ValidationError`).
- `prepare.cli.base.read_accessions_file` becomes a thin re-export (keeps the CLI import path).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_accessions.py
import pytest
from prepare.core.exceptions import ValidationError
from prepare.data.accessions import read_accessions_file, require_accession


def test_reader_skips_comments_blanks_and_duplicates(tmp_path):
    f = tmp_path / "acc.txt"
    f.write_text("# header\nSRR1\n\n  SRR2  \nSRR1\n")
    assert read_accessions_file(f) == ["SRR1", "SRR2"]


@pytest.mark.parametrize("bad", ["../../victim", "..", "/tmp/x", "SRR1/..", "SRR 1", "srr1"])
def test_reader_refuses_unsafe_lines(tmp_path, bad):
    f = tmp_path / "acc.txt"
    f.write_text(f"SRR1\n{bad}\n")
    with pytest.raises(ValidationError) as info:
        read_accessions_file(f)
    assert "line 2" in str(info.value) and bad in str(info.value)


def test_require_accession():
    assert require_accession("ERR123") == "ERR123"
    with pytest.raises(ValidationError):
        require_accession("../x")
```

```python
# tests/test_data_sra.py (new test in the download_sra class)
def test_download_sra_refuses_traversal_before_touching_disk(self, tmp_path):
    victim = tmp_path / "victim"
    victim.mkdir()
    (tmp_path / "acc.txt").write_text("../../victim\n")
    fastq = tmp_path / "proj" / "fastq"
    with pytest.raises(ValidationError):
        download_sra(fastq, tmp_path / "acc.txt")
    assert victim.is_dir()
    assert not (tmp_path / "proj").exists()
```

- [ ] **Step 2: Run them**

Run: `conda run -n metaquest python -m pytest tests/test_accessions.py tests/test_data_sra.py -k "traversal or accessions" -v`
Expected: FAIL (ImportError, then the victim folder is removed).

- [ ] **Step 3: Implement**

```python
# prepare/data/accessions.py
"""The accessions file: one SRA run accession per line."""
import re
from pathlib import Path
from typing import List, Union

from prepare.core.constants import SRA_ACCESSION_PATTERN
from prepare.core.exceptions import DataAccessError, ValidationError

_ACCESSION_RE = re.compile(SRA_ACCESSION_PATTERN)


def require_accession(accession: str) -> str:
    """``accession`` when it is a well-formed SRA run accession; ``ValidationError`` otherwise.

    The pattern admits only upper-case letters and digits, so a value that passes is safe
    to use as a folder name.
    """
    if not _ACCESSION_RE.fullmatch(accession or ""):
        raise ValidationError(f"Invalid SRA accession format: {accession!r}")
    return accession


def read_accessions_file(path: Union[str, Path], validate: bool = True) -> List[str]:
    """Accessions listed in ``path`` in file order, without blanks, comments or repeats."""
    try:
        with open(path, "r") as handle:
            lines = handle.read().splitlines()
    except OSError as e:
        raise DataAccessError(f"Cannot read accessions file {path}: {e}") from e
    seen: dict = {}
    for number, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if validate and not _ACCESSION_RE.fullmatch(line):
            raise ValidationError(f"{path}, line {number}: {line!r} is not an SRA accession")
        seen.setdefault(line, None)
    return list(seen)
```

In `prepare/cli/base.py` replace the body of `read_accessions_file` with `from prepare.data.accessions import read_accessions_file as _read; return _read(path)` (keep the name exported; `accessions_from_args` unchanged). In `download.py:323-324` replace the open loop with `all_accessions = read_accessions_file(accessions_file)`. In `metadata.py` make `_read_accessions_file` call the shared reader. At the top of `download_accession` (after the docstring) add `require_accession(accession)`; the same first line in `link_dataset` (`prepare/store/link.py:201`) and `sra_dir` (`prepare/store/layout.py:122`), which are the choke points for Task 1.2.

- [ ] **Step 4: Run the tests and the suite slices**

Run: `conda run -n metaquest python -m pytest tests/test_accessions.py tests/test_data_sra.py tests/test_cli_commands.py tests/test_store_link.py tests/test_store_layout.py -q`
Expected: PASS. If a fixture uses a non-conforming accession such as `"acc1"` in a store path, change the fixture to `SRR1`-style names (record the files touched in the commit message).

- [ ] **Step 5: Commit**

```bash
git add prepare/data/accessions.py prepare/cli/base.py prepare/data/sra/download.py prepare/data/sra/accession.py prepare/data/metadata.py prepare/store/link.py prepare/store/layout.py tests/test_accessions.py tests/test_data_sra.py
git commit -m "fix: validate accessions before they become paths (audit H1)"
```

### Task 1.2: Validate accessions in store_link, store_unlink, store_verify (M1)

**Files:**
- Modify: `prepare/cli/commands/store/link.py:139-160, 212-237`, `prepare/cli/commands/store/verify.py:85-88`
- Test: `tests/test_cli_store.py`

**Interfaces:** Consumes `require_accession` from Task 1.1. The CLI commands validate every positional accession up front and return exit code 2 with one error line listing the bad values, before any filesystem or registry action.

- [ ] **Step 1: Failing tests** (copy the `_link_args` / `_unlink_args` / `_verify_args` helper pattern at `tests/test_cli_store.py:91-119`)

```python
class TestStoreCommandsRefuseUnsafeAccessions:
    def test_unlink_outside_fastq_is_refused(self, tmp_path, monkeypatch, caplog):
        root = tmp_path / "store"; init_store(root)
        outside = tmp_path / "important_link"; outside.symlink_to(tmp_path)
        proj = tmp_path / "proj"; (proj / "fastq").mkdir(parents=True)
        monkeypatch.chdir(proj)
        rc = StoreUnlinkCommand().execute(_unlink_args(["../important_link"], data_root=str(root)))
        assert rc == 2 and outside.is_symlink()
        assert "not an SRA accession" in caplog.text

    def test_link_outside_sra_is_refused(self, tmp_path, monkeypatch):
        root = tmp_path / "store"; init_store(root)
        proj = tmp_path / "proj"; proj.mkdir(); monkeypatch.chdir(proj)
        rc = StoreLinkCommand().execute(_link_args(["../tmp"], data_root=str(root), accept_partial=True))
        assert rc == 2 and not (proj / "fastq").exists()
        assert not (proj / "prepare_registry.json").exists()

    def test_verify_refuses_bad_accession(self, tmp_path):
        root = tmp_path / "store"; init_store(root)
        assert StoreVerifyCommand().execute(_verify_args(["../x"], data_root=str(root))) == 2
```

- [ ] **Step 2: Run** `pytest tests/test_cli_store.py -k Unsafe -v`; expected FAIL.

- [ ] **Step 3: Implement.** Add to `prepare/cli/base.py` a helper the three commands share:

```python
def refuse_bad_accessions(logger: logging.Logger, accessions: Sequence[str]) -> Optional[int]:
    """Exit code 2 after one error line when any accession is malformed; ``None`` otherwise."""
    bad = [acc for acc in accessions if not re.fullmatch(SRA_ACCESSION_PATTERN, acc)]
    if bad:
        logger.error("not an SRA accession: %s", ", ".join(bad))
        return 2
    return None
```

Call it first thing in `StoreLinkCommand.execute`, `StoreUnlinkCommand.execute` and `StoreVerifyCommand.execute` (`if (rc := refuse_bad_accessions(self.logger, args.accessions)) is not None: return rc`). Keep the `require_accession` guards added in Task 1.1 in `link_dataset` and `sra_dir` as defence in depth.

- [ ] **Step 4: Run** `pytest tests/test_cli_store.py -q`; expected PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: store link, unlink and verify refuse malformed accessions (audit M1)"`

### Task 1.3: Resolve the registry path before it reaches the catalogue; live links protect datasets (H2)

**Files:**
- Modify: `prepare/data/registry.py:162-163`, `prepare/store/usage.py:30-32, 181-201`, `prepare/cli/commands/store/gc.py:243-257`
- Test: `tests/test_store_usage.py`, `tests/test_store_gc.py`, `tests/test_data_registry.py`

**Interfaces:**
- `registry_path(explicit)` returns `Path(explicit).resolve()`.
- `_upsert_project_row` stores `str(Path(registry.path).resolve())`.
- `linked_by(paths, catalog, accession) -> List[str]` no longer skips stale projects: a symlink on this host that resolves to the dataset protects it whatever the registry state. The name is suffixed `" (stale)"` when the project is stale, so gc's reason line says so.

- [ ] **Step 1: Failing tests**

```python
# tests/test_data_registry.py, in TestLoadSave
def test_registry_path_resolves_an_explicit_relative_path(self, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert registry_path("prepare_registry.json") == (tmp_path / "prepare_registry.json").resolve()

# tests/test_store_usage.py
def test_relative_registry_is_stored_absolute(tmp_path, paths, monkeypatch):
    monkeypatch.chdir(tmp_path)
    registry = _registry(tmp_path, "p1", "proj", str(tmp_path))
    registry.path = Path("prepare_registry.json")
    with catalog_write(paths) as cat:
        _upsert_project_row(cat, registry, "p1")
        row = cat.conn.execute("SELECT registry FROM projects").fetchone()
    assert row["registry"] == str((tmp_path / "prepare_registry.json").resolve())

# tests/test_store_usage.py, TestLinkedBy: replace the "stale project's symlink is ignored" test with
def test_a_stale_projects_live_symlink_still_counts(self, ...):
    ... same setup as the old test ...
    assert linked_by(paths, cat, "SRR1") == ["proj (stale)"]

# tests/test_store_gc.py (end-to-end, from the audit reproduction)
def test_gc_include_stale_keeps_a_dataset_linked_by_a_relative_registry_project(tmp_path, monkeypatch, capsys):
    root = tmp_path / "store"; init_store(root)
    acc_dir = _write_dataset_dir(store_paths(root), "SRR2")
    ... write sidecar + upsert_dataset as in tests/test_store_gc.py:378-396 ...
    proj = tmp_path / "proj"; proj.mkdir(); monkeypatch.chdir(proj)
    assert StoreInitCommand().execute(_init_args(root, proj, registry="prepare_registry.json")) == 0
    assert StoreLinkCommand().execute(_link_args(["SRR2"], data_root=str(root), registry="prepare_registry.json")) == 0
    used = root / "locks" / "SRR2.used"; os.utime(used, (1, 1))  # past the 24 h grace
    monkeypatch.chdir(tmp_path)
    rc = StoreGcCommand().execute(_gc_args(data_root=str(root), include_stale=True, yes=True, json=True))
    assert rc == 0 and acc_dir.is_dir()
    report = json.loads(capsys.readouterr().out)
    assert "SRR2" not in [c["accession"] for c in report.get("removed", [])]
```

- [ ] **Step 2: Run** the four tests; expected FAIL (the gc test removes `sra/SRR2`).
- [ ] **Step 3: Implement** the three one-line changes plus the `linked_by` loop: build `stale_names = {row["project_id"]: row for row in stale_projects(catalog)}` once, do not `continue` on stale ids, and append `f"{name} (stale)"` when the id is stale. In `_classify_dataset` the "linked" branch already reads `linked_by`, so no change there.
- [ ] **Step 4: Run** `pytest tests/test_store_usage.py tests/test_store_gc.py tests/test_data_registry.py tests/test_cli_store.py -q`; expected PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: resolve registry paths in the catalogue; live links protect datasets from gc (audit H2)"`

### Task 1.4: Project folders come from the registry location, not the cwd (H4 and the store_link/unlink/identity variants)

**Files:**
- Modify: `prepare/data/registry.py:53-61` (`ProjectPaths`), `prepare/cli/commands/status/command.py:170-178`, `prepare/cli/commands/store/link.py:139-160, 212-237`, `prepare/store/usage.py:49-59`
- Test: `tests/test_cli_status.py`, `tests/test_cli_store.py`, `tests/test_data_registry.py`

**Interfaces:**
- `ProjectPaths.anchored(cls, root: Path, fastq="fastq", metadata="metadata", genomes="genomes", targeted="targeted", matches="matches") -> ProjectPaths`: each folder that is relative is joined onto `root`; absolute folders pass through.
- `project_root(registry_file: Path) -> Path` returns `registry_file.resolve().parent`.
- `ensure_project_identity` records `project.path` as `project_root(registry.path)` when the registry has a path, else cwd.

- [ ] **Step 1: Failing tests**

```python
# tests/test_cli_status.py
def test_reconcile_from_a_subfolder_uses_the_project_above(tmp_path, monkeypatch):
    root = _project_tree(tmp_path)           # helper at line 240 writes fastq/SRR000001 and the registry
    sub = root / "sub"; sub.mkdir(); monkeypatch.chdir(sub)
    assert StatusCommand().execute(_args(registry=None, reconcile=True)) == 0
    reg = load_registry(root / "prepare_registry.json")
    assert reg.datasets["SRR000001"]["download"]["state"] == "downloaded"
    assert not (sub / "fastq").exists()

# tests/test_cli_store.py
def test_store_link_from_a_subfolder_links_into_the_project_fastq(...):
    ... store with SRR1 as at tests/test_cli_store.py:1545-1572; registry at proj/; chdir proj/sub ...
    assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root))) == 0
    assert (proj / "fastq" / "SRR1").is_symlink() and not (proj / "sub" / "fastq").exists()

def test_ensure_project_identity_records_the_registry_folder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path / "elsewhere")  # create it first
    registry = Registry(path=tmp_path / "proj" / "prepare_registry.json")
    assert ensure_project_identity(registry)["path"] == str((tmp_path / "proj").resolve())
```

- [ ] **Step 2: Run**; expected FAIL.
- [ ] **Step 3: Implement.**

```python
# prepare/data/registry.py
def project_root(registry_file: Union[str, Path]) -> Path:
    """The project folder: where its registry file lives."""
    return Path(registry_file).resolve().parent


@dataclass
class ProjectPaths:
    ...
    @classmethod
    def anchored(cls, root: Path, **folders: Union[str, Path]) -> "ProjectPaths":
        """Folders joined onto ``root`` unless already absolute."""
        values = {}
        for name in ("fastq", "metadata", "genomes", "targeted", "matches"):
            raw = Path(folders.get(name, name))
            values[name] = raw if raw.is_absolute() else root / raw
        return cls(**values)
```

In status `execute`: compute `registry_file = registry_path(args.registry)` first, then `paths = ProjectPaths.anchored(project_root(registry_file), fastq=args.fastq_folder, ...)`. In store link and unlink: `fastq = Path(args.fastq_folder); if not fastq.is_absolute(): fastq = project_root(registry_path(args.registry)) / fastq` and use `fastq` everywhere `args.fastq_folder` was used (including `record_download(..., fastq_dir)`). In `ensure_project_identity`: `base = project_root(registry.path) if registry.path else Path.cwd().resolve()`; name and path from `base`.

- [ ] **Step 4: Run** `pytest tests/test_cli_status.py tests/test_cli_store.py tests/test_data_registry.py tests/test_registry_reconcile.py tests/test_registry_reconcile_reverify.py -q`; expected PASS. Existing tests pass absolute folders, so they are unaffected.
- [ ] **Step 5: Commit** `git commit -m "fix: status, store_link and store_unlink act on the project that owns the registry (audit H4)"`

### Task 1.5: Registry forward compatibility and malformed JSON (likely: registry)

**Files:** `prepare/data/registry.py:176-231`; test `tests/test_data_registry.py`.

**Interfaces:** `load_registry` raises `DataAccessError` for a non-object document, a non-integer version, or a `version` above `SCHEMA_VERSION` ("written by a newer prepare; upgrade prepare or use that version"). Unknown top-level keys are kept in a new `Registry.extra: Dict[str, Any]` and written back.

- [ ] **Step 1: Tests**: `[]` document raises `DataAccessError`; `{"version": 3}` raises `DataAccessError` mentioning "newer"; `{"version": "x"}` raises `DataAccessError`; a registry with `"notes": {"a": 1}` round-trips through `save_registry` with the key intact.
- [ ] **Step 2: Run**; FAIL with AttributeError / TypeError / key lost.
- [ ] **Step 3: Implement**: after `json.loads`, `if not isinstance(data, dict): raise DataAccessError(...)`; parse version with `try: version = int(data.get("version", SCHEMA_VERSION)) except (TypeError, ValueError)`; `if version > SCHEMA_VERSION: raise DataAccessError(...)`; add `extra: Dict[str, Any] = field(default_factory=dict)` holding every key outside the seven known; `_write_registry` builds `payload = {**registry.extra, "version": ..., ...}`.
- [ ] **Step 4: Run** `pytest tests/test_data_registry.py tests/test_cli_status.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: the registry refuses newer schemas and keeps unknown keys"`

### Task 1.6: Link takes the dataset lock; gc copes with a missing folder (likely: link/gc race, catalogue drift)

**Files:** `prepare/store/link.py:201-254`, `prepare/cli/commands/store/gc.py:506-531`, `prepare/data/registry_reconcile.py:240-270`; tests `tests/test_store_link.py`, `tests/test_store_gc.py`, `tests/test_registry_reconcile.py`.

**Interfaces:**
- `link_dataset(..., locked: bool = False)`: unless `locked` (the store hand-off already holds the lock at `store_handoff.py:521`), wraps the `is_dir` check, symlink or copytree in `dataset_lock(paths, accession, blocking=False)`; `LockHeld` becomes `DataAccessError("<acc> is being written by <holder>; retry later")`.
- gc `_remove_one_dataset`: `FileNotFoundError` on `os.replace` of the dataset folder deletes the catalogue row inside the same `catalog_write` pass (add the accession to a `vanished` list that `_remove_candidates` deletes without an aside) and reports it under `removed` with reason "folder already gone".
- `_mark_missing` calls `update_linked(registry, acc, add=False)` for each accession it marks missing.

- [ ] **Step 1: Tests**: (a) with `dataset_lock(paths, "SRR1")` held in the test, `link_dataset(..)` raises `DataAccessError` and leaves no link; (b) the store hand-off path still links while holding the lock (`locked=True`); (c) a catalogue row for `SRR9` with no folder: `store_gc --yes` removes the row and the report says "folder already gone"; (d) reconcile on a registry with `store.linked == ["SRR1"]` and no `fastq/SRR1` leaves `store.linked` empty.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement** as described; `store_handoff.py:521` passes `locked=True`.
- [ ] **Step 4: Run** `pytest tests/test_store_link.py tests/test_store_gc.py tests/test_registry_reconcile.py tests/test_data_sra.py -k "store or link or reconcile" -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: store links take the dataset lock; gc drops rows whose folder is gone; reconcile clears store.linked"`

### Task 1.7: Wave 1 gate

- [ ] Run `conda run -n metaquest make check` and `conda run -n metaquest make test`; both must pass. Fix anything that fails inside this task and commit as `test: wave 1 gate fixes`.

---

## Wave 2: wrong results

### Task 2.1: PERMANOVA (H3)

**Files:** `prepare/processing/diversity.py:312-363`; tests `tests/test_diversity.py`.

**Interfaces:** `_permanova_test(distance_matrix, groups, n_permutations, rng=None) -> Tuple[float, float]` implementing Anderson (2001): drop samples whose group label is NaN (and the matching rows/columns); `SS_T = sum_{i<j} d_ij^2 / N`; `SS_W = sum_g sum_{i<j in g} d_ij^2 / n_g`; `SS_A = SS_T - SS_W`; `F = (SS_A/(a-1)) / (SS_W/(N-a))`; `p = (1 + #{F_perm >= F}) / (1 + n_permutations)`; `F = nan` and `p = nan` when `a < 2` or `N - a < 1`. `rng` is a `numpy.random.Generator` so tests are deterministic.

- [ ] **Step 1: Failing tests**

```python
def _reference_f(d, labels):
    """Plain-loop PERMANOVA pseudo-F, Anderson 2001 eq. 4, written independently of the code."""
    n = len(labels); groups = sorted(set(labels)); a = len(groups)
    ss_t = sum(d[i][j] ** 2 for i in range(n) for j in range(i + 1, n)) / n
    ss_w = 0.0
    for g in groups:
        idx = [i for i in range(n) if labels[i] == g]
        ss_w += sum(d[i][j] ** 2 for i in idx for j in idx if i < j) / len(idx)
    return ((ss_t - ss_w) / (a - 1)) / (ss_w / (n - a))


def test_permanova_matches_the_reference_formula():
    rng = np.random.default_rng(0)
    points = np.vstack([rng.normal(0, 1, (5, 3)), rng.normal(2, 1, (5, 3))])
    d = squareform(pdist(points))
    labels = ["a"] * 5 + ["b"] * 5
    f, p = _permanova_test(d, pd.Series(labels), 199, rng=np.random.default_rng(1))
    assert f == pytest.approx(_reference_f(d, labels))
    assert 0 < p <= 1 and p >= 1 / 200


def test_permanova_two_identical_clusters_has_f_near_one_and_high_p(): ...  # permuted labels on one cluster, p > 0.2

def test_permanova_nan_labels_are_dropped_and_single_group_is_nan():
    d = np.ones((4, 4)) - np.eye(4)
    f, p = _permanova_test(d, pd.Series(["a", "a", np.nan, np.nan]), 10)
    assert np.isnan(f) and np.isnan(p)
```

- [ ] **Step 2: Run** `pytest tests/test_diversity.py -k permanova -v`; FAIL (F off by the subtracted means).
- [ ] **Step 3: Implement** the formula above with `np.triu_indices`; permutation loop uses `rng.permutation(labels)`; `perform_permanova` passes `rng=np.random.default_rng(seed)` where `seed` is a new optional argument defaulting to `None`. Keep the output keys unchanged.
- [ ] **Step 4: Run** `pytest tests/test_diversity.py tests/test_cli_commands_advanced_analysis.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: PERMANOVA uses the standard sums of squares and the +1 p-value (audit H3)"`

### Task 2.2: Threshold grid and header validation in Branchwater parsing (M6, M7)

**Files:** `prepare/data/branchwater.py:134-148, 344-347, 372-381`; tests `tests/test_data_branchwater.py`.

- [ ] **Step 1: Tests**: `_threshold_counts(np.array([0.3, 0.6, 0.7]), 0.1)` gives counts `[0,0,0,1,2,3,3,3,3,3,3]` for thresholds `1.0 ... 0.0`; `process_branchwater_files` on a folder holding `bad.csv` with header `foo,bar` copies nothing, appends `"bad.csv"` to `errors`, and the existing `test_process_branchwater_files_invalid_format` is rewritten to use a real file instead of mocking `validate_header`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: `raw = np.round(np.arange(int(round(1 / step_size)), -1, -1) * step_size, 10)` and compare `ordered` against `np.round(ordered, 10)`; at lines 136 and 347 call `BranchWaterFormatPlugin._validate_headers(csv_file, headers)` (which raises `ValidationError`, a `PrepareError`, already caught by both `except` clauses). Rename `_validate_headers` to a public `require_headers` and keep `validate_header` for callers that want a bool.
- [ ] **Step 4: Run** `pytest tests/test_data_branchwater.py tests/test_containment_pin.py tests/test_plugins.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: containment thresholds are exact decimals; bad Branchwater headers are refused (audit M6, M7)"`

### Task 2.3: NCBI taxonomy rate limit, cache dtypes, BOM (M8, L2, L4)

**Files:** `prepare/data/taxonomy.py:84, 328-331`, `prepare/plugins/formats/branchwater.py:101`; tests `tests/test_taxonomy_extended.py:152-159`, `tests/test_taxonomy.py`, `tests/test_plugins.py`.

- [ ] **Step 1: Tests**: `NCBITaxonomyClient("a@b", api_key=None).request_delay == pytest.approx(1/3)` and with a key `== pytest.approx(1/10)` (replace `test_api_key_faster_rate_limit`'s assertion direction); a cache file whose first row has an empty `tax_id` followed by a row `562`: `validate_taxonomic_assignments(["Escherichia coli"], ..., cache_file=...)` returns `tax_id == "562"` and name columns are strings not NaN; `BranchWaterFormatPlugin.parse_file` on a file starting with `﻿acc,containment\n` parses one row.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: swap the delay expression; `pd.read_csv(cache_file, encoding="utf-8", dtype=str, keep_default_na=False)`; `open(file_path, "r", encoding="utf-8-sig", newline="")`.
- [ ] **Step 4: Run** `pytest tests/test_taxonomy.py tests/test_taxonomy_extended.py tests/test_plugins.py tests/test_plugins_comprehensive.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: NCBI pacing without a key is 3 requests per second; taxonomy cache keeps strings; BOM tolerated (audit M8, L2, L4)"`

### Task 2.4: Single-run metadata files and unchecked NCBI bodies (M10, likely: 200 error body)

**Files:** `prepare/data/metadata.py:198-228`; tests `tests/test_data_metadata.py`.

**Interfaces:** `_download_single_metadata` runs the response through `_split_efetch_packages(xml, {accession})`; if the accession is absent, returns `(False, "NCBI returned no package for <acc>")` and writes nothing; otherwise writes the single-run package.

- [ ] **Step 1: Tests**: patch `Entrez.efetch` to return a two-run package (SRR1 then SRR2) when asked for SRR2: the written `SRR2_metadata.xml` parses with `_extract_metadata_fields` to `Run_ID == "SRR2"`; patch it to return `<html>error</html>`: no file is written and the result is `(False, ...)`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement** as described; non-XML bodies raise `ET.ParseError` inside `_split_efetch_packages`, catch and report.
- [ ] **Step 4: Run** `pytest tests/test_data_metadata.py tests/test_cli_metadata.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: single-run metadata downloads keep only the requested run (audit M10)"`

### Task 2.5: Forced re-extraction that maps nothing removes the old reads (M9)

**Files:** `prepare/data/extraction_staging.py` (new `remove_published(out_dir, genome_id) -> List[Path]`), `prepare/data/read_extraction.py:461-480`; tests `tests/test_extraction_staging.py`, `tests/test_read_extraction.py`.

- [ ] **Step 1: Tests**: `remove_published` deletes the files named in `.<G>.extracted.json` and the record, returns the removed paths, and is a no-op without a record; in `_map_and_extract` (patch `_run_minimap2`/`_count_records` to give zero), pre-existing `targeted/<ACC>/<G>_1.fastq.gz` and the record are gone afterwards.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: both zero branches call `remove_published(out_dir, genome_id)` before returning.
- [ ] **Step 4: Run** `pytest tests/test_extraction_staging.py tests/test_read_extraction.py tests/test_cli_read_extraction.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: a forced re-extraction with no mapped reads removes the earlier FASTQs (audit M9)"`

### Task 2.6: Quality profiling on zero-length reads and 3' homopolymers (M5, L1)

**Files:** `prepare/sra/quality.py:181, 280-287`; tests `tests/test_sra_analytics.py`.

- [ ] **Step 1: Tests**: `_sample_head` on a FASTQ holding one empty record and one `ACGT` record returns `n_contents == [0.0, 0.0]`; `_analyze_sequence_complexity(["ACGTAAAA"])["homopolymer_rate_mean"] == pytest.approx(1 / 8)`; `["AAAACGTA"]` gives the same.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: `n_contents.append(sequence.count("N") / len(sequence) if sequence else 0.0)`; after the loop `if current_run >= 4: homopolymer_count += 1`.
- [ ] **Step 4: Run** `pytest tests/test_sra_analytics.py tests/test_cli_sra_profile.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: head sampler tolerates empty reads; trailing homopolymer runs count (audit M5, L1)"`

### Task 2.7: Pass --min-abundance through (L3)

**Files:** `prepare/data/taxonomy.py:495-520`, `prepare/cli/commands/advanced_analysis.py:394-410`; tests `tests/test_taxonomy.py`, `tests/test_cli_commands_advanced_analysis.py`.

- [ ] Add `min_abundance: float = 0.001` to `analyze_taxonomic_composition`, forward it to `create_taxonomic_summary`; the command passes `min_abundance=args.min_abundance`. Test: with `min_abundance=0.5` a 0.1-abundance taxon is absent from the genus summary; with the default it is present. Commit `git commit -m "fix: taxonomic_summary honours --min-abundance (audit L3)"`.

### Task 2.8: Likely data-layer bugs with cheap fixes

**Files and changes, each with one test in the module's existing test file:**
- `prepare/data/gtdb.py:94-97, 183-190`: `search_taxon` pages (`page` parameter, loop until a short page, cap 50 pages); test with a patched `requests.get` returning 100 then 3 rows gives 103 accessions.
- `prepare/data/genome_taxonomy.py:115-120`: use `rows[0]` only if `rows[0].get("accession") == accession` (also accept the `GCA_`/`GCF_` sibling); else `None`. Test: a row for a different accession yields `None`.
- `prepare/data/sra_metadata.py:161-169`: search term built as `" OR ".join(f"{acc}[ACCN]" for acc in accessions)` and `retmax = str(len(accessions) * 2)`; after the fetch, log the accessions absent from the result at WARNING. Test: the term contains `[ACCN]` and the missing accession is logged.
- `prepare/data/branchwater_search.py:267-269`: do not write the cache when `len(kept) == 0` unless the response is non-empty (`text.strip()` has rows); log "empty result not cached". Test: a zero-row response leaves no cache file.
- `prepare/data/genome_download.py:144-184`: `extract_dir = zip_path.parent / f"ncbi_extract_{os.getpid()}"`; sort the glob; after the loop compare requested accessions (new parameter `expected: Optional[Sequence[str]]`) with the extracted set and raise `DataAccessError` listing the missing ones. Test: a zip with one of two requested accessions raises and names the missing one.
- `prepare/data/defaults.py:49-55`: strip `.gz` before looking at the suffix (`Path(path.name.removesuffix(".gz")).suffix`). Test: `_separator_for(Path("a.tsv.gz")) == "\t"`.
- `prepare/processing/selection.py:77-79, 93`: `sort_values(ascending=False, kind="stable")` after a `sort_index()` so ties resolve by name. Test: two equal values return in index order.
- `prepare/data/metadata.py:562-564`: pick the `SRAFile` whose `semantic_name == "SRA Normalized"` (fallback: the last one). Test: the normalized URL comes from the matching element regardless of order.
- `prepare/data/taxonomy.py:371-398`: an empty species list returns an empty frame with the `_VALIDATION_CACHE_COLUMNS` and logs 0 of 0, with no division. Test: `validate_taxonomic_assignments([], "a@b")` returns an empty DataFrame.
- `prepare/processing/diversity.py:167-168`: samples missing from `metadata` raise `ValidationError` naming them before any file is written. Test: missing sample gives `ValidationError`.

- [ ] One commit per bullet or grouped by file, messages prefixed `fix:` with the audit's wording.

### Task 2.9: Wave 2 gate

- [ ] `make check` and `make test`; fix and commit `test: wave 2 gate fixes`.

---

## Wave 3: download result handling

### Task 3.1: Carry `force` through the retry pass (M2)

**Files:** `prepare/data/sra/retry.py:169-258`, `prepare/data/sra/download.py` (the call site of `_retry_failed_downloads`); tests `tests/test_data_sra.py:1713-1731`.

**Interfaces:** `_retry_failed_downloads(..., force: bool = False, ...)` passes `force=force` to the downloader. `download_sra` passes its own `force`.

- [ ] **Step 1: Tests**: replace `test_retry_failed_downloads_uses_force_false` with `test_retry_passes_the_callers_force` asserting `call_args.kwargs["force"] is True` when called with `force=True`, and `False` by default. Add an end-to-end test: `download_sra(..., force=True, max_retries=1)` with a downloader that fails once then succeeds never reports "already exists" (patch `download_accession`).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `pytest tests/test_data_sra.py tests/test_stop_token.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: the retry pass keeps --force (audit M2)"`

### Task 3.2: One place strips the retry prefix (M3)

**Files:** `prepare/data/sra/accession.py` (move `_RETRY_PREFIX_RE` and `_body` here as public `result_body(message) -> str`), `prepare/data/sra/run_report.py:52-61` (import it), `prepare/data/sra/retry.py:157-162`, `prepare/cli/commands/sra.py:515-531`; tests `tests/test_sra_run_report.py`, `tests/test_cli_commands.py`, `tests/test_data_sra.py`.

- [ ] **Step 1: Tests**: `_split_not_found(["SRR1"], {"SRR1": "Retry 1: incomplete: 3 of 9 spots"})` puts SRR1 in `not_found`; the result recorder given `("SRR1", True, "Retry 1: linked from store")` records `source == "store"`, `attempt` not counted, `store.linked` contains SRR1 and a usage row is queued; given `"Retry 1: already exists"` it takes the `_record_present` path.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: in `_record_result` set `body = result_body(message)` and use `body` in the three comparisons (keep `message` for the registry text); in `_split_not_found` use `result_body(str(message))`.
- [ ] **Step 4: Run** `pytest tests/test_sra_run_report.py tests/test_cli_commands.py tests/test_data_sra.py tests/test_cli_sra_verdicts.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: retry-prefixed results are classified like first-pass results (audit M3)"`

### Task 3.3: Space guard and duplicate accessions (M4)

**Files:** `prepare/data/sra/space.py:209-258`; tests `tests/test_sra_space.py`.

- [ ] Reader de-duplication from Task 1.1 removes the realistic trigger; make the guard itself safe: `reserve` raises `RuntimeError(f"{accession} is already reserved")` if `accession in self._held` (a programming error, not a user error); test that two reserves without a release raise and that reserve, release, reserve, release leaves `_reserved` at zero. Commit `git commit -m "fix: the space guard refuses a double reservation (audit M4)"`.

### Task 3.4: Stop token inside the retry loop (likely)

**Files:** `prepare/data/sra/retry.py:242-276`; tests `tests/test_stop_token.py`.

- [ ] At the top of the per-accession loop: `if accession_mod.stop_requested(stop): failed_accessions.extend(retry_batch[index:]); break` without writing "interrupted" results or notifying. Test: a stop set after the first accession leaves the second's `download_results` untouched and `on_result` called once. Commit `git commit -m "fix: a stop during the retry pass leaves the remaining accessions unrecorded"`.

### Task 3.5: Wave 3 gate

- [ ] `make check`, `make test`; commit fixes.

---

## Wave 4: reports and plots

### Task 4.1: Metadata colouring and sizing in PCA and t-SNE (H5)

**Files:** `prepare/visualization/interactive.py:48-75, 233-240, 404-412`; tests `tests/test_interactive_plots.py`.

- [ ] **Step 1: Tests**: with a metadata frame indexed `["s2", "s1", "s3"]` and `grp = ["b", "a", "a"]`, `create_interactive_pca(data, metadata, color_by="grp")` returns a figure whose trace names are `{"a", "b"}` (not `"nan"`); `size_by` numeric column is applied (marker sizes vary); the same for `create_interactive_tsne`; `create_diversity_comparison_plot(alpha, None, "grp")` raises `ValidationError("--metadata-file is required ...")`; a sample absent from the metadata raises `ValidationError` naming it.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: `plot_data[color_by] = metadata.loc[sample_names, color_by].to_numpy()` (same for size and t-SNE); add a shared `_metadata_values(metadata, sample_names, column) -> np.ndarray` that checks `column in metadata.columns`, finds `missing = set(sample_names) - set(metadata.index)` and raises `ValidationError` naming the missing samples; use it in the diversity plot with an explicit `None` check first.
- [ ] **Step 4: Run** `pytest tests/test_interactive_plots.py tests/test_interactive_plots_extended.py tests/test_cli_commands_advanced_analysis.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: PCA and t-SNE colour and size by metadata again (audit H5)"`

### Task 4.2: Explorer escaping and metadata columns (M11, M12)

**Files:** `prepare/visualization/explorer.py:118-121, 297-298, 361`; tests `tests/test_explorer.py`.

- [ ] **Step 1: Tests**: a taxonomy family `"<img src=x onerror=alert(2)>"` appears escaped (`&lt;img`) in the output and never raw; metadata with columns `geo_loc_name_country_calc` and `collection_date_sam` fills the Location and Date cells; a NaN value renders as an empty cell, not `"nan"`; a duplicated metadata index does not raise.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement**: `jinja2.Environment(loader=jinja2.BaseLoader(), autoescape=True)` and remove the manual `html_escape(title)` (the `|safe` filters on pre-built HTML stay); define `_LOCATION_COLUMNS = ("geo_loc_name_country_calc", "geo_loc_name_sam", "location")` and `_DATE_COLUMNS = ("collection_date_sam", "collection_date")`, take the first present; `metadata = metadata[~metadata.index.duplicated(keep="first")]`; `"" if pd.isna(value) else str(value)`.
- [ ] **Step 4: Run** `pytest tests/test_explorer.py tests/test_cli_explore.py -q`; PASS.
- [ ] **Step 5: Commit** `git commit -m "fix: the explorer escapes family names and reads the metadata columns prepare writes (audit M11, M12)"`

### Task 4.3: Explorer scaling (gap)

**Files:** `prepare/visualization/explorer.py:96-126, 232-255`; tests `tests/test_explorer.py`.

- [ ] Vectorise `_build_long_dataframe` with `containment_df.stack()` plus a merge on sample id; cap the table at `EXPLORER_TABLE_ROWS = 5000` (new constant in `prepare/core/constants.py`) sorted by containment, with a note in the header "showing top N of M rows"; heatmap height capped at 2000 px. Test: 300 samples by 50 genomes builds in under 2 s and the table holds at most 5000 rows. Commit `git commit -m "perf: the explorer builds its table with pandas and caps the rows shown"`.

### Task 4.4: Map coordinate parsing and the unused country column (M13, gap)

**Files:** `prepare/plugins/visualizers/map.py:67-122, 216-285`; tests `tests/test_plugins_comprehensive.py`.

- [ ] Parser: regex `^\s*(-?\d+(?:\.\d+)?)\s*([NS])?\s*[, ]\s*(-?\d+(?:\.\d+)?)\s*([EW])?\s*$`; non-string input returns `(None, None)`; range check (|lat| <= 90, |lon| <= 180) else `(None, None)`. `country_column`: when both coordinates are missing and the column is present, log once how many rows had only a country and skip them (no geocoding). Tests: `"38.98 N 77.10 W"` gives `(38.98, -77.10)`; `"52.62 N 1.22 E"`, `"12.5,-3.1"`, `"-33.9 S 151.2 E"`, `3.5` (float) and `"999 N 0 E"` give the expected values or `(None, None)`. Commit `git commit -m "fix: the map plugin parses NCBI lat_lon strings (audit M13)"`.

### Task 4.5: sourmash plugin errors and partial parses (M14, M15)

**Files:** `prepare/plugins/sourmash_plugin.py:25-75, 115-132, 174-225, 264-301`; tests `tests/test_sourmash_plugin.py`.

- [ ] **Interfaces:** each plugin `main` first calls `setup_logging(level=logging.INFO)` from `prepare.utils.logging` (console handler on stderr); `prepare_parse` gains `--details-file`, `--registry`, `--registry-max-screened` and mirrors `cli/commands/containment.py:70-93` (collect `errors`, record screening in the registry, exit 1 with the error line when any file failed). Tests: running `PrepareParsePlugin.main` with an unreadable match file exits 1 and the message appears in `capsys.readouterr().err`; `PreparePlotPlugin.main` with a missing file prints "Error plotting containment" to stderr.
- [ ] Commit `git commit -m "fix: the sourmash plugins report errors and record screening like the CLI (audit M14, M15)"`

### Task 4.6: Static plots: pie colours, save location, title, figure closing (L5, L6, likely: --title, plt.close)

**Files:** `prepare/visualization/plots.py:151-156, 197-200, 277, 338-373`, `prepare/cli/commands/containment.py:128`, `prepare/cli/commands/metadata.py:366-391`, `prepare/plugins/visualizers/{bar,heatmap,map}.py` savefig sites; tests `tests/test_visualization_plots.py`, `tests/test_cli_commands.py`.

- [ ] Changes: `colors = [colors] if isinstance(colors, str) else colors` at the top of `plot_metadata_counts` and `plot_containment`; `plot_metadata_counts` saves through `plot_output_path`-style paths beside the input (`Path(file_path).parent / f"{stem}_{plot_type}.{fmt}"`) and skips saving for DataFrame input; the radar fallback passes `limit`; `--title` given without `--show-title` shows the title (`show_title = args.show_title or args.title is not None` in both commands); every `savefig` is followed by `plt.close(fig)` in a `finally`, and `_save_plot_if_needed` returns the path written. Without `--save-format` the metadata command logs at WARNING "no --save-format given; nothing was written". Tests: `--colors red` on a pie draws; `counts_bar.png` lands in the input's folder; the title text appears when only `title=` is passed; `len(plt.get_fignums()) == 0` after `plot_containment(..., save_format="png")`.
- [ ] Commit `git commit -m "fix: metadata plots save beside their input, accept one colour, keep a given title and close figures (audit L5, L6)"`

### Task 4.7: Report HTML: stat card, CSV export, sort order (gap, likely)

**Files:** `prepare/visualization/project_report.py:124`, `prepare/utils/html.py:310-359`; tests `tests/test_project_report.py`, `tests/test_html_assets.py`.

- [ ] Add `"analysed"` to the stat-card loop; `exportCSV` prefixes a cell starting with `=`, `+`, `-`, `@` with a single quote; the comparator parses a number only when the whole trimmed cell is numeric (`/^-?\d+(\.\d+)?$/`), numbers sort before text, and the first click sorts ascending. Tests: the report HTML contains `>analysed<`; the script text contains the formula guard regex and `/^-?\d+/`. Commit `git commit -m "fix: project report shows the analysed stage; CSV export guards formulas; table sort is total (audit gaps)"`.

### Task 4.8: Wave 4 gate

- [ ] `make check`, `make test`; commit fixes.

---

## Wave 5: CLI polish, CI parity, packaging

### Task 5.1: Broken pipes, index column, -v clash, dry run (M16, L7, L8, L9)

**Files:** `prepare/cli/main.py:245-258`, `prepare/cli/commands/explore.py:254-299`, `prepare/cli/base.py:138-146`, `prepare/cli/commands/store/status.py:47`, `prepare/data/metadata.py:284-307`; tests `tests/test_exit_codes.py`, `tests/test_cli_explore.py`, `tests/test_logging_setup.py`, `tests/test_data_metadata.py`.

- [ ] Changes: in `_run_command` add `except BrokenPipeError: devnull = os.open(os.devnull, os.O_WRONLY); os.dup2(devnull, sys.stdout.fileno()); return 0` before the generic handler; `find_by_taxonomy` writes with `index=False` (and the usage error case returns 2 via `parser.error` semantics: raise `argparse.ArgumentError` is not available in `execute`, so return the usage exit code constant from `prepare.core.exceptions`; grep for it and add `EXIT_USAGE = 2` beside the other codes if none exists); `store_status` renames its own flag to `--details` (keep `--verbose` as a hidden alias that logs a deprecation line) so the global `-v` is always wired; `download_metadata` checks `dry_run` before `ensure_directory`. Tests: a `BrokenPipeError` raised by a command gives exit 0; the TSV's first header cell is `Run_ID` (or the index name); `main(["store_status", "-v", ...])` parses; `download_metadata(..., dry_run=True)` creates no folder.
- [ ] Commit `git commit -m "fix: broken pipes exit 0, find_by_taxonomy drops the index column, -v works everywhere, dry runs write nothing (audit M16, L7, L8, L9)"`

### Task 5.2: Settings and logging edge cases (likely: memory grammar, --log-level case, flag conflicts, doctor fallback, JSON startup errors, store config types, nan floats)

**Files:** `prepare/core/settings.py:200-208`, `prepare/cli/base.py:107-112`, `prepare/cli/main.py:195-203, 297-313`, `prepare/store/resolve.py:91-105`; tests `tests/test_core_settings.py`, `tests/test_logging_setup.py`, `tests/test_exit_codes.py`, `tests/test_store_resolve.py`.

- [ ] Changes: `_memory` delegates to `resources.parse_memory(value, limit=None)` inside `try/except ValueError` so startup and assembly agree; `--log-level` gets `type=str.upper`; `--log-level X` together with `-v`/`-q` is a usage error (`parser.error`) like `-q -v`; the doctor fallback copies `log_file`, `log_quiet`, `log_verbose` from `parsed_args` onto the default runtime; when `--json` is set and startup fails, write `{"error": "..."}` to stdout before returning 3; `resolve_store_root` validates `store` is a table and `data_root` a string (`ConfigurationError` otherwise), reads the config only after the explicit and env candidates miss, and `expanduser()`s candidates; float settings reject `nan`/`inf` at parse time with a `ValueError` (usage code 2 through argparse `type=`). One test per bullet.
- [ ] Commit `git commit -m "fix: settings grammar, log level handling and store config types (audit likely bugs)"`

### Task 5.3: CI parity and release checks (gaps, L10, L11)

**Files:** `.github/workflows/ci.yml`, `.github/workflows/release.yml`, `Makefile:58`, `pyproject.toml`, `prepare/py.typed` (new, empty), `scripts/check_version.py` (new), `scripts/check_docs_commands.py:58-72`, `scripts/check_no_print.sh`; tests `tests/test_check_docs_commands.py` (new), `tests/test_packaging.py` (new).

- [ ] Changes:
  - `ci.yml` lint job: replace the black, flake8 and mypy steps with `make check` (it already includes them) after installing `radon`; the test job unchanged; the build job adds `python scripts/check_version.py` and `unzip -l dist/*.whl | grep prepare/py.typed`.
  - `release.yml` build job: add a step `python scripts/check_version.py --tag "${GITHUB_REF_NAME}"` before `python -m build`, and a `pytest tests/ -q -x` step after installing `.[dev]`.
  - `scripts/check_version.py`: reads `pyproject.toml` `version` and `prepare/__init__.py` `__version__`, exits 1 if they differ or (with `--tag vX.Y.Z`) differ from the tag.
  - `Makefile:58`: name `plotly_js_script()`.
  - `check_docs_commands.py`: a command counts as documented only when it appears as a table row ``| `name` |`` or an invocation `prepare name` (also indented or after `$ `); also scan `docs/*.md` for unknown invocations.
  - `check_no_print.sh`: pattern also matches `sys\.stdout\.buffer\.write` and `from sys import stdout`; exclude `print(..., file=sys.stderr)` is not needed because no such call exists (verify with grep before committing).
  - `tests/test_packaging.py`: `prepare/py.typed` exists; versions agree.
- [ ] Commit `git commit -m "ci: run make check in CI, check the release tag against the version, ship py.typed (audit gaps, L10, L11)"`

### Task 5.4: Help text and cache names (gaps)

- [ ] Store-true options stop printing "(default: False)" (use `argparse.SUPPRESS` default rendering or strip in the formatter in `prepare/cli/base.py`); `--verify-downloads` keeps one default note; `taxonomy_validation` writes `taxonomy_validation_cache.csv` and `enrich_taxonomy` keeps `taxonomy_cache.tsv` with the docs naming both (`docs/analysis_commands.md`). Test: `prepare download_sra --help` output has no "(default: False)". Commit `git commit -m "cli: cleaner help defaults and distinct taxonomy cache names"`.

### Task 5.5: Wave 5 gate

- [ ] `make check`, `make test`, `make build`; commit fixes.

---

## Wave 6: structural

### Task 6.1: One NCBI request pacer

**Files:** create `prepare/data/ncbi_pacing.py`; modify `prepare/data/taxonomy.py:72-119`, `prepare/data/sra_metadata.py:59-112`, `prepare/data/metadata.py:41-70`; tests `tests/test_ncbi_pacing.py` (new) plus the three clients' tests.

**Interfaces:** `class RequestPacer: def __init__(self, api_key: Optional[str], per_second_without_key=3, per_second_with_key=10)`, `def wait(self) -> None` (thread-safe, monotonic clock), `def delay(self) -> float`. The three sites construct one pacer each and call `wait()` before each request; the module-level `_last_request_time`, `_PACE_LOCK` and `request_delay` attributes are removed (tests that read `request_delay` move to `pacer.delay`).

- [ ] Tests: `RequestPacer(None).delay == pytest.approx(1/3)`; with a key `1/10`; two `wait()` calls at least `delay` apart (patch `time.monotonic`/`time.sleep`). Commit `git commit -m "refactor: one NCBI request pacer for the three clients"`.

### Task 6.2: Finish the plugin base

**Files:** `prepare/plugins/base.py`, `prepare/plugins/__init__.py`, `prepare/visualization/plots.py:20-28`, `prepare/data/branchwater.py:30-35`, `prepare/core/constants.py:281`; tests `tests/test_plugins.py`.

- [ ] Changes: `FormatPlugin(Plugin)` with abstract `parse_file(cls, path, genome_id)` and `validate_header`; `VisualizerPlugin(Plugin)` with abstract `create_plot(self, data, **kwargs)`; `PluginRegistry.register` is idempotent for the same class object and raises only for a different class under the same name; registration moves from import side effects into `prepare.plugins.load_builtin_plugins()` called once from `prepare/plugins/__init__.py` plus `discover_entry_points(group=PLUGIN_REGISTRY_NAME)` that loads `importlib.metadata.entry_points(group=...)`, logging and skipping a plugin that fails to import; delete the leftover TypeVars and comment at lines 125-131; drop the `# type: ignore[attr-defined]` at the call sites. Tests: re-registering the same class is a no-op; a different class with the same name raises; a fake entry point that raises on load is skipped with a warning; `importlib.reload(prepare.visualization.plots)` does not raise.
- [ ] Commit `git commit -m "refactor: abstract plugin interfaces, idempotent registration and entry-point discovery"`

### Task 6.3: A map command, or drop the maps extra

- [ ] Add `prepare plot_map --metadata-file <tsv> --output <png> [--lat-lon-column lat_lon]` in `prepare/cli/commands/metadata.py` (group: metadata) using `MapVisualizerPlugin`; it `require()`s the `maps` extra; add the command to the README table and `docs/analysis_commands.md` (the docs-command gate enforces this; command count becomes 45). Test: with cartopy patched out, the command returns the `MissingExtra` exit code with the install hint; with a stub `create_plot` it writes the file. Commit `git commit -m "feat: plot_map draws sample locations from metadata, so the maps extra has a command"`.

### Task 6.4: Shared helpers for duplicated logic

- [ ] `prepare/utils/sizes.py`: `tree_bytes(path) -> int` (follows no symlinks, `is_file(follow_symlinks=False)` via `os.lstat`, OSError-safe) replaces gc `_path_bytes`, link `_tree_bytes` and the inner loop of `cleanup.transient_bytes`; `prepare/visualization/interactive.py`: `_prepare_matrix(data, scale) -> Tuple[np.ndarray, List[str]]` shared by PCA and t-SNE and `_write_figure(fig, output_file)` used by the four `write_html` sites, writing through `atomic_path`; `prepare/data/gtdb.py`: one `_write_cache` used by both writers; `prepare/cli/commands/genome.py`: use `read_accessions_file(validate=False)` for genome accession files (GCF/GCA pattern is different; add `pattern=` parameter to the reader). Tests: `tree_bytes` ignores a symlink that points outside; `_write_figure` leaves no temp file on failure. Commit `git commit -m "refactor: shared size walker, matrix preparation and figure writer"`.

### Task 6.5: gc performance and journal pruning

- [ ] `linked_by(paths, catalog, accession, stale_ids=None)`: callers pass a precomputed stale set; `_dataset_candidates` computes `stale_projects` once per run and `_still_removable` once per removal batch; `store_gc` prunes `usage.jsonl` lines for removed accessions and `store_reindex` compacts `projects.jsonl` to the latest line per project. Catalogue: `migrate()` refuses a `schema_version` above its own with `DataAccessError`. Tests: `stale_projects` is called once for a 50-dataset classification (count via monkeypatch); after gc the journal has no lines for the removed accession; a catalogue with `schema_version = 99` is refused. Commit `git commit -m "perf: gc loads project registries once; journals are pruned; newer catalogues refused"`.

### Task 6.6: Module splits

- [ ] Split without behaviour change, each with its re-exports so imports in tests keep working for one release:
  - `prepare/data/registry.py` (1001): move `scan_extractions`, `bootstrap_from_disk` and the disk-scanning helpers (lines ~800-1001) to `prepare/data/registry_scan.py`; move `load_registry`, `_write_registry`, `save_registry`, `registry_transaction`, `registry_update` to `prepare/data/registry_io.py`; `registry.py` keeps the dataclasses, `registry_path`, `project_root`, record functions and re-exports.
  - `prepare/data/read_extraction.py` (960): `_export_mapped_fastq`, `reference_coverage`, `coverage_table_path` and the samtools/seqkit helpers move to `prepare/data/read_extraction_export.py`.
  - `prepare/cli/commands/sra.py` (775): `_result_recorder`, `_record_present`, `_sidecar_completeness` move to `prepare/cli/commands/sra_record.py` as module functions taking the command's logger.
  - `prepare/data/metadata.py` (797): the Entrez download functions (`_download_single_metadata`, `_download_batch_metadata`, `_split_efetch_packages`, `_download_accessions_individually`, pacing) move to `prepare/data/metadata_fetch.py`.
  - `tests/test_data_sra.py` (4358): split by class into `test_data_sra_download.py`, `test_data_sra_retry.py`, `test_data_sra_store.py` (mechanical, keep fixtures in a `tests/_sra_fixtures.py`).
  - Add to `tests/test_module_sizes.py` a frozen ceiling of 700 lines for the split modules so they do not regrow.
- [ ] Gate: `make check` (module size check) and `make test`. One commit per module: `refactor: split <module> by responsibility`.

### Task 6.7: Rename the CSS prefix and the remaining leftovers

- [ ] Replace `mq-` with `pre-` in `prepare/utils/html.py`, `prepare/visualization/explorer.py`, `prepare/sra/reporting.py`, `prepare/visualization/project_report.py`, `tests/test_html_assets.py:21`, `tests/test_explorer.py:63` (`sed -i '' 's/mq-/pre-/g'` on those six files, then grep that no `mq-` remains under `prepare/` or `tests/`); historical plans under `docs/superpowers/plans/` are left as they are (history). Test: `test_html_assets` asserts the new class names. Commit `git commit -m "style: report CSS classes use the pre- prefix"`.

### Task 6.8: Release 0.10.1

- [ ] `prepare/__init__.py` and `pyproject.toml` to `0.10.1`; CHANGELOG `## [0.10.1] - <date>` with sections Fixed (every H/M/L bug in one line each, user-facing wording), Changed (accession validation at the CLI, stale projects' links now protect datasets, PERMANOVA p-values change, `store_status --details`, `plot_map` added, `pre-` CSS prefix), Internal (pacer, plugin base, splits). Update `docs/ARCHITECTURE.md` module list for the new files and `docs/analysis_commands.md` for `plot_map` and PERMANOVA notes. Run `make check`, `make test`, `make build`, and the fresh-venv wheel install smoke used for 0.10.0. Commit `git commit -m "chore: release 0.10.1"`.
- [ ] Hand-off: `superpowers:finishing-a-development-branch` (merge to main, tag `v0.10.1`, push; the Release and CI workflows must be green).

---

## Verification

End to end, after Task 6.8 on a scratch project:
1. `printf 'SRR1\n../../x\n' > acc.txt && prepare download_sra --accessions-file acc.txt --dry-run` exits 2 with a line naming line 2, and nothing is created.
2. Store flow from the audit reproduction: `store_init --registry prepare_registry.json`, `store_link SRR2`, age the `.used` marker, `store_gc --include-stale --yes` from another folder: the dataset stays and the report says "still linked by proj (stale)".
3. `cd proj/sub && prepare status --reconcile`: the parent registry keeps its downloaded states and no `sub/fastq` appears.
4. `prepare diversity_analysis ... --permanova-formula` on `test_data/` gives an F matching the test's reference function for the same input.
5. `prepare interactive_plot --plot-type pca --color-by <col>` on `test_data/` produces traces named after the column's values.
6. `prepare find_by_taxonomy --family X | head -1` exits 0 and the first header cell is not empty.
7. `sourmash scripts prepare_plot --file-path /nonexistent` prints an error line to stderr.
8. `make check`, `make test` (3300+ tests), `make build` green locally; CI, Release, Branchwater and smoke workflows green on `v0.10.1`.
