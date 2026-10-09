# Rename to PREPARE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the tool from MetaQuest to PREPARE: package, command, distribution, on-disk names, repository, folder and documentation, released as 0.10.0.

**Architecture:** One scripted, mechanical rename commit (Task 1) makes the tree import under the new name; three implementers then clean up disjoint areas in parallel (Tasks 2-4), a reviewer checks the whole branch (Task 5), the gates run (Task 6), and the outward-facing steps follow after the merge (Task 7).

**Tech Stack:** Python 3.12, git, GitHub CLI, conda env `metaquest` (holds the editable install; it keeps its name).

**Spec:** `docs/superpowers/specs/2026-10-09-rename-to-prepare-design.md`

## Global Constraints

- Package, command, PyPI distribution and conda environment name in the Makefile are all `prepare`.
- Prose uses `PREPARE`; identifiers use `Prepare` (`PrepareError`, `PrepareParsePlugin`).
- Settings variables are `PREPARE_<NAME>`; config path is `~/.config/prepare/config.toml`.
- On-disk names: `prepare_registry.json`, `prepare_store.json`, `prepare_version`, `.prepare/runs/`, `.prepare-assembly.json`, `fastq/.prepare-tmp/`, thread name `prepare-lock-heartbeat`, `PLUGIN_REGISTRY_NAME = "prepare_plugins"`, sourmash entry points `prepare_parse`, `prepare_plot`, `prepare_diversity`, `prepare_taxonomy`.
- No fallback to any old name; no compatibility shim.
- Version 0.10.0 in `pyproject.toml` and `prepare/__init__.py`.
- Untouched: `.git/`, `.claude/`, `.remember/`, `.superpowers/`, `docs/superpowers/` (every file, including this plan and its spec), `CHANGELOG.md` below the 0.10.0 entry, the old command names in `prepare/cli/commands/renamed.py`.
- ASCII only in `prepare/`, `tests/`, `scripts/`, `Makefile`, `setup.cfg`, `pyproject.toml`; no Unicode anywhere new.
- Modest scientific language in documentation and code.
- Every commit message ends with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

## Review Focus

1. A project folder from 0.9.0 (holding `metaquest_registry.json`) run under `prepare status`: expected behaviour is "no registry here" with the usual hint, not a crash and not silent reading of the old file. Test pinned to Task 2.
2. `METAQUEST_DATA` still set in the user's shell after the upgrade: `prepare` must ignore it and find no store (clean break), and `prepare doctor` must not mention it. Test pinned to Task 2.
3. The `metaquest` console script left in the environment by the old editable install: after `pip install -e .` the `prepare` script exists; the stale `metaquest` script fails to import and is removed by the reinstall. Covered by Task 3 step on reinstall and by Task 6.
4. The sourmash plugin entry points: `sourmash scripts prepare_parse` must resolve to `PrepareParsePlugin`; a stale `metaquest_parse` entry must not remain in `pyproject.toml`. Test pinned to Task 2.
5. The word "prepare" already occurs in prose and code before the rename (for example "prepare the folder"); the mechanical substitution must not produce doubled or garbled text such as "prepareprepare". Residual search in Task 5 covers `metaquest`; Task 5 also greps for `prepareprepare` and `PREPAREPREPARE`.

---

### Task 1: Mechanical rename commit (main session)

**Files:**
- Delete worktrees: `.claude/worktrees/reporting`, `.claude/worktrees/storage-fixes`
- Rename: `metaquest/` to `prepare/`
- Modify: every tracked text file outside the untouched set
- Create (scratchpad only, not committed): `rename_prepare.py`

**Interfaces:**
- Produces: the branch `rename/prepare` whose head imports as `prepare`, with the substitutions below applied uniformly. Tasks 2-4 start from this commit.

- [ ] **Step 1: Remove the stale worktrees and create the branch**

```bash
cd /Users/andreassjodin/Code/metaquest
git worktree remove .claude/worktrees/reporting
git worktree remove .claude/worktrees/storage-fixes
git branch -D feat/reporting fix/storage-resume
git worktree prune
git checkout -b rename/prepare
```

Expected: `git worktree list` shows only the main checkout.

- [ ] **Step 2: Write the rename script in the scratchpad**

```python
"""Mechanical rename of MetaQuest to PREPARE over the tracked files of a git checkout."""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve()
SKIP_PREFIXES = ("docs/superpowers/", ".claude/", ".remember/", ".superpowers/", ".git/")
SKIP_FILES = {"CHANGELOG.md"}
MARKDOWN = {".md"}

# Longer and upper-case forms first. In prose a bare MetaQuest/Metaquest becomes PREPARE;
# one followed by an identifier character (MetaQuestError, MetaquestParsePlugin) becomes Prepare.
RULES = [
    (re.compile(r"METAQUEST"), "PREPARE"),
    (re.compile(r"Meta[Qq]uest(?=[A-Za-z0-9_])"), "Prepare"),
    (re.compile(r"Meta[Qq]uest"), "PREPARE"),
    (re.compile(r"metaquest"), "prepare"),
]


def tracked() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [line for line in out.splitlines() if line]


def wanted(rel: str) -> bool:
    if rel in SKIP_FILES or rel.startswith(SKIP_PREFIXES):
        return False
    return True


def rewrite(text: str) -> str:
    for pattern, replacement in RULES:
        text = pattern.sub(replacement, text)
    return text


def main() -> None:
    changed = 0
    for rel in tracked():
        if not wanted(rel):
            continue
        path = ROOT / rel
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        new = rewrite(text)
        if new != text:
            path.write_text(new, encoding="utf-8")
            changed += 1
    print(f"rewrote {changed} files")


if __name__ == "__main__":
    main()
```

Save as `$SCRATCHPAD/rename_prepare.py`.

- [ ] **Step 3: Move the package and run the script**

```bash
cd /Users/andreassjodin/Code/metaquest
git mv metaquest prepare
python -I "$SCRATCHPAD/rename_prepare.py" .
python -I "$SCRATCHPAD/rename_prepare.py" . >/dev/null   # second run must change 0 files
```

Expected: the first run reports about 300 files rewritten; the second reports `rewrote 0 files`.

- [ ] **Step 4: Rename AGENTS.md content in place without adding it**

```bash
python - <<'PY'
import re, pathlib
p = pathlib.Path("AGENTS.md"); t = p.read_text()
for pat, rep in [(r"METAQUEST","PREPARE"),(r"Meta[Qq]uest(?=[A-Za-z0-9_])","Prepare"),(r"Meta[Qq]uest","PREPARE"),(r"metaquest","prepare")]:
    t = re.sub(pat, rep, t)
p.write_text(t)
PY
git status --short | grep AGENTS.md
```

Expected: `?? AGENTS.md` (still untracked).

- [ ] **Step 5: Confirm the tree imports and the fast tests run**

```bash
conda run -n metaquest pip install -e . -q
conda run -n metaquest prepare --help | head -3
conda run -n metaquest python -m pytest tests/test_optional_imports.py tests/test_core_settings.py -q -x
```

Expected: the help text's first line names `prepare`; both test files pass.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor: mechanical rename of the metaquest package and names to prepare

git mv metaquest prepare plus a scripted substitution: METAQUEST to PREPARE,
MetaQuest/Metaquest to PREPARE in prose and Prepare in identifiers, metaquest to
prepare. CHANGELOG.md and docs/superpowers/ are untouched.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: Package code and name constants (Opus)

**Files:**
- Modify: everything under `prepare/` that the script produced; especially `prepare/core/constants.py`, `prepare/core/settings.py`, `prepare/core/exceptions.py`, `prepare/store/resolve.py`, `prepare/store/sidecar.py`, `prepare/store/catalog.py`, `prepare/data/registry.py`, `prepare/data/run_log.py`, `prepare/data/assembly_identity.py`, `prepare/utils/lockfile.py`, `prepare/plugins/sourmash_plugin.py`, `prepare/cli/main.py` (parser description), `prepare/cli/commands/doctor.py`
- Test: `tests/test_core_settings.py`, `tests/test_registry*.py`, `tests/test_store*.py`, `tests/test_run_log.py`, `tests/test_assembly_*.py`, and a new module `tests/test_rename_clean_break.py`

**Interfaces:**
- Consumes: the Task 1 commit.
- Produces: `prepare.data.registry.REGISTRY_FILENAME == "prepare_registry.json"`, `prepare.core.constants.STORE_MARKER == "prepare_store.json"`, `PLUGIN_REGISTRY_NAME == "prepare_plugins"`, settings prefix `PREPARE_`, and the class names `PrepareError`, `PrepareParsePlugin`, `PreparePlotPlugin`, `PrepareDiversityPlugin`, `PrepareTaxonomyPlugin`. Task 3 relies on these class names in `pyproject.toml`.

- [ ] **Step 1: Read the spec's name-mapping table and grep the package for each new name**

```bash
grep -rn "prepare_registry.json\|prepare_store.json\|prepare_version\|\.prepare/\|\.prepare-assembly\|\.prepare-tmp\|prepare-lock-heartbeat\|prepare_plugins\|PREPARE_" prepare | wc -l
grep -rn "PrepareError\|PrepareParsePlugin\|PreparePlotPlugin\|PrepareDiversityPlugin\|PrepareTaxonomyPlugin" prepare | wc -l
```

Expected: both counts are above zero; read each hit once and fix any grammar the substitution broke (for example a docstring "PREPARE's" is fine, "the PREPARE package" is fine, "Prepare" standing alone in prose is wrong and becomes "PREPARE").

- [ ] **Step 2: Write the failing clean-break tests**

Create `tests/test_rename_clean_break.py`:

```python
"""The 0.10.0 rename is a clean break: old MetaQuest names are not read."""
from importlib.metadata import entry_points

from prepare.core import constants
from prepare.core.settings import SETTINGS
from prepare.data.registry import REGISTRY_FILENAME, load_registry
from prepare.store.resolve import resolve_store_root


def test_registry_filename_is_prepare():
    assert REGISTRY_FILENAME == "prepare_registry.json"


def test_store_marker_is_prepare():
    assert constants.STORE_MARKER == "prepare_store.json"


def test_plugin_registry_name_is_prepare():
    assert constants.PLUGIN_REGISTRY_NAME == "prepare_plugins"


def test_every_setting_env_var_uses_prepare_prefix():
    for spec in SETTINGS.values():
        for name in (spec.env, *spec.aliases):
            assert name.startswith("PREPARE_"), name


def test_old_registry_file_is_not_read(tmp_path):
    (tmp_path / "metaquest_registry.json").write_text('{"datasets": {"SRR1": {}}}')
    registry = load_registry(tmp_path)
    assert registry.path == tmp_path / "prepare_registry.json"
    assert not (tmp_path / "prepare_registry.json").exists()
    assert registry.datasets == {}  # use the real attribute that holds the dataset map


def test_old_store_variable_is_ignored(tmp_path, monkeypatch):
    monkeypatch.setenv("METAQUEST_DATA", str(tmp_path / "old"))
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    assert resolve_store_root(explicit=None, registry_root=None) is None


def test_sourmash_entry_points_use_prepare_names():
    names = {ep.name for ep in entry_points(group="sourmash.cli_script")}
    assert "prepare_parse" in names
    assert not any(name.startswith("metaquest") for name in names)
```

The comment on the dataset-map assertion means: open `prepare/data/registry.py`, use the `Registry` attribute that holds the per-accession entries (it may be named differently), and delete the comment. `load_registry` returns an empty registry bound to the new file name when that file is absent; the test checks it never reads the old one.

- [ ] **Step 3: Run the tests; fix the names so they pass**

```bash
conda run -n metaquest python -m pytest tests/test_rename_clean_break.py -v
```

Expected after fixes: 7 passed. The entry-point test needs `pip install -e .` to have run after Task 1 (it did, in Task 1 step 5); if it fails only because the metadata is stale, run `conda run -n metaquest pip install -e . -q` once and rerun.

- [ ] **Step 4: Run the owned test modules**

```bash
conda run -n metaquest python -m pytest tests/test_core_settings.py tests/test_registry_assembly.py tests/test_registry_blocks.py tests/test_registry_concurrent_commands.py tests/test_registry_reconcile.py tests/test_registry_reconcile_reverify.py tests/test_registry_timing.py tests/test_run_log.py tests/test_assembly_identity.py tests/test_assembly_layout.py tests/test_store_adopt.py tests/test_store_catalog.py tests/test_store_gc.py tests/test_store_journal.py tests/test_store_layout.py tests/test_store_link.py tests/test_store_locks.py tests/test_store_resolve.py tests/test_store_sidecar.py tests/test_store_stats.py tests/test_store_usage.py tests/test_store_verdicts.py -q
```

Expected: all pass.

- [ ] **Step 5: Lint and type-check the package**

```bash
conda run -n metaquest python -m black prepare tests/test_rename_clean_break.py
conda run -n metaquest python -m flake8 prepare tests/test_rename_clean_break.py
conda run -n metaquest python -m mypy prepare
```

Expected: no output from flake8, mypy reports no issues.

- [ ] **Step 6: Commit**

```bash
git add prepare tests/test_rename_clean_break.py tests/test_core_settings.py tests/test_registry*.py tests/test_store*.py tests/test_run_log.py tests/test_assembly_*.py
git commit -m "refactor: package code uses the PREPARE names; clean-break tests

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: Tests, scripts, packaging and CI (Sonnet)

**Files:**
- Modify: every file under `tests/` not owned by Task 2, `scripts/check_no_print.sh`, `scripts/check_ascii.sh`, `scripts/check_docs_commands.py`, `scripts/check_atomic_writes.sh`, `scripts/atomic_writes_allowlist.txt`, `scripts/smoke_chain.sh`, `Makefile`, `setup.cfg`, `pyproject.toml` (all lines except `version`), `environment.yml`, `requirements*.txt`, `local_test.sh`, `.github/workflows/ci.yml`, `.github/workflows/release.yml`, `.github/workflows/audit.yml`, `.github/workflows/branchwater.yml`, `.github/workflows/smoke.yml`, `tests/test_module_sizes.py`, `tests/helpers_processes.py`, `tests/test_performance_regressions.py`

**Interfaces:**
- Consumes: the class names from Task 2 (`PrepareParsePlugin` and friends) in the `pyproject.toml` entry points; the Task 1 commit.
- Produces: the quality gates (`make check`) pointing at `prepare/`; CI workflows named `PREPARE CI` and `PREPARE Release`.

- [ ] **Step 1: Check the packaging metadata**

```bash
grep -n "prepare" pyproject.toml
```

Expected lines include `name = "prepare"`, `prepare = "prepare.cli.main:main"`, the four `prepare_*` sourmash entry points pointing at `prepare.plugins.sourmash_plugin:Prepare*Plugin`, `include = ["prepare*"]`, `source = ["prepare"]`, the mypy override `module = ["prepare.data.registry_blocks"]`, and the perf marker text naming `PREPARE_PERF_SCALE`. The URLs point at `github.com/FOI-Bioinformatics/prepare`. Fix any line the script missed.

- [ ] **Step 2: Check the scripts and the Makefile**

```bash
grep -n "prepare" scripts/check_no_print.sh scripts/check_ascii.sh scripts/check_docs_commands.py scripts/check_atomic_writes.sh Makefile setup.cfg environment.yml
```

Expected: the paths read `prepare/`, the conda environment in the Makefile is `prepare`, `environment.yml` has `name: prepare`. Fix the help text grammar of the Makefile if the substitution produced "the 'prepare' conda environment" (that is correct; leave it).

- [ ] **Step 3: Check the CI workflows**

```bash
grep -n "prepare\|PREPARE" .github/workflows/*.yml
```

Expected: `name: PREPARE CI`, `name: PREPARE Release`; the pip-audit comment reads `unpublished "prepare" package`; every `black`, `flake8`, `mypy`, `pytest --cov` target is `prepare`. Fix anything missed.

- [ ] **Step 4: Check the process helpers and perf tests**

```bash
grep -n "prepare\|PREPARE" tests/helpers_processes.py tests/test_performance_regressions.py tests/test_module_sizes.py | head -40
```

Expected: `python -m prepare.cli.main`, `PREPARE_DATA` in the isolated environment, `PREPARE_PERF_SCALE`, and `KNOWN_EXCEPTIONS` keys `prepare/data/read_extraction.py` and `prepare/data/registry.py`.

- [ ] **Step 5: Reinstall and run the gate scripts**

```bash
conda run -n metaquest pip install -e ".[dev]" -q
conda run -n metaquest bash scripts/check_no_print.sh
conda run -n metaquest bash scripts/check_ascii.sh
conda run -n metaquest bash scripts/check_atomic_writes.sh
conda run -n metaquest python scripts/check_docs_commands.py
conda run -n metaquest python tests/test_module_sizes.py --check
```

Expected: each script exits 0. `check_docs_commands.py` may fail until Task 4 finishes the README; if its only complaint is a README command, note it in the report and move on.

- [ ] **Step 6: Run the tests you own**

```bash
conda run -n metaquest python -m pytest tests/ -q -x --deselect tests/test_rename_clean_break.py -p no:cacheprovider
```

Expected: all pass (Task 2's modules are also run; a failure there is reported, not fixed here).

- [ ] **Step 7: Commit**

```bash
git add tests scripts Makefile setup.cfg pyproject.toml environment.yml requirements*.txt local_test.sh .github
git reset -q pyproject.toml && git add -p pyproject.toml   # stage everything but the version line, or leave the version line unchanged (Task 4 bumps it)
git commit -m "build: tests, scripts, packaging and CI use the prepare name

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

If `git add -p` is impractical, leave the `version` line at 0.9.0 and stage the whole file; Task 4 changes only that line.

---

### Task 4: Documentation, CHANGELOG and version (Sonnet)

**Files:**
- Modify: `README.md`, `docs/README.md`, `docs/ARCHITECTURE.md`, `docs/branchwater_workflow.md`, `docs/configuration.md`, `docs/hpc.md`, `docs/packaging.md`, `docs/pipeline_overview.md`, `docs/SRA_ENHANCED_FEATURES.md`, `CONTRIBUTING.md`, `CLAUDE.md`, `CHANGELOG.md` (intro line and a new 0.10.0 entry only), `pyproject.toml` (the `version` line only), `prepare/__init__.py` (`__version__` only)

**Interfaces:**
- Consumes: the Task 1 commit; the on-disk names from the Global Constraints.
- Produces: version `0.10.0`; a CHANGELOG entry that the release notes use.

- [ ] **Step 1: Read every documentation file once for grammar the substitution broke**

Look for "Prepare" standing alone in prose (wrong, use "PREPARE"), "PREPARE's" (fine), headings, badge URLs, `pip install prepare`, `prepare --help`, `sourmash scripts prepare_parse`, `~/.config/prepare/config.toml`, `PREPARE_<NAME>` and the file names in the Global Constraints. The README's first paragraph describes PREPARE as a command-line toolkit for analyzing metagenomic datasets based on genome containment; PREPARE has no expansion. Do not touch `docs/superpowers/`.

- [ ] **Step 2: Add the CHANGELOG entry**

Change the intro line to `All notable changes to PREPARE are documented in this file.` and insert above `## [0.9.0]`:

```markdown
## [0.10.0] - 2026-10-09

### Changed

- The tool is renamed from MetaQuest to PREPARE, because the old name clashes with an unrelated
  consumer product. The package, the command, the distribution and the repository
  (`github.com/FOI-Bioinformatics/prepare`) are now `prepare`; the exception base class is
  `PrepareError`; the sourmash plugin entry points are `prepare_parse`, `prepare_plot`,
  `prepare_diversity` and `prepare_taxonomy`.
- Every name the tool reads or writes follows the rename, with no fallback to the old name. A
  project created with 0.9.0 or earlier is used with 0.10.0 by renaming, by hand:
  `metaquest_registry.json` to `prepare_registry.json`; `.metaquest/` to `.prepare/`;
  each `.metaquest-assembly.json` marker to `.prepare-assembly.json`; `fastq/.metaquest-tmp/`
  to `fastq/.prepare-tmp/` (or delete it); in a shared store, `metaquest_store.json` to
  `prepare_store.json` in the store root and in every dataset folder; the `metaquest_version`
  key inside the registry and the sidecars to `prepare_version`; `~/.config/metaquest/config.toml`
  to `~/.config/prepare/config.toml`; and every `METAQUEST_<NAME>` environment variable to
  `PREPARE_<NAME>` (`METAQUEST_DATA` to `PREPARE_DATA`, `METAQUEST_PERF_SCALE` to
  `PREPARE_PERF_SCALE`). After renaming the store files, run `prepare store_reindex`.
- The hidden pointers for the commands renamed in 0.5.0 (`sra_stats`, `sra_profile_quality`,
  `sra_dashboard`, `sra_compare`) are unchanged.
```

Before writing this entry, confirm each file name against the code (`grep -rn "prepare_store.json\|prepare_version" prepare/store prepare/data`) and against the spec; if the sidecar and catalogue names differ in the code, state each real name.

- [ ] **Step 3: Bump the version**

```bash
sed -i '' 's/^version = "0.9.0"/version = "0.10.0"/' pyproject.toml
sed -i '' 's/^__version__ = "0.9.0"/__version__ = "0.10.0"/' prepare/__init__.py
grep -n '0.10.0' pyproject.toml prepare/__init__.py
```

Expected: one hit in each file. Update the version named in CLAUDE.md's architecture paragraph ("As of 0.9.0 ... 44 commands") to 0.10.0 only if the sentence states the current version; leave historical references.

- [ ] **Step 4: Run the documentation gate**

```bash
conda run -n metaquest python scripts/check_docs_commands.py
grep -rn "MetaQuest\|metaquest\|METAQUEST" README.md docs/*.md CONTRIBUTING.md CLAUDE.md
```

Expected: the gate passes; the grep prints nothing.

- [ ] **Step 5: Commit**

```bash
git add README.md docs/*.md CONTRIBUTING.md CLAUDE.md CHANGELOG.md pyproject.toml prepare/__init__.py
git commit -m "docs: PREPARE documentation, 0.10.0 changelog entry and version bump

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: Whole-branch review and residual search (Opus)

**Files:**
- Read: the whole diff `git diff de4a953..HEAD`
- Modify: none (findings go to the main session, which assigns fixes to the owning task's implementer)

- [ ] **Step 1: Residual search**

```bash
grep -rIn -i "metaquest" . --exclude-dir=.git --exclude-dir=.claude --exclude-dir=.remember --exclude-dir=.superpowers --exclude-dir=superpowers --exclude=CHANGELOG.md
grep -rIn "prepareprepare\|PREPAREPREPARE\|PreparePrepare" . --exclude-dir=.git --exclude-dir=.claude
awk '/^## \[0.9.0\]/{exit} {print}' CHANGELOG.md | grep -in "metaquest"
```

Expected: the first grep prints only the clean-break tests (`tests/test_rename_clean_break.py`, which name the old file and variable on purpose) and nothing else; the second prints nothing; the third prints only lines inside the 0.10.0 entry.

- [ ] **Step 2: Review the diff against the spec**

Check: every name in the spec's mapping table appears in the new form; no behaviour change beyond names; prose reads as PREPARE, identifiers as Prepare; the ASCII rule holds; commit messages carry the trailer. Report findings as Critical / Important / Minor with file and line.

---

### Task 6: Quality gates and build (main session)

- [ ] **Step 1: Gates**

```bash
cd /Users/andreassjodin/Code/metaquest
conda run -n metaquest make check
conda run -n metaquest make test
conda run -n metaquest make pipeline
conda run -n metaquest make build
```

Expected: each exits 0; `dist/` holds `prepare-0.10.0-py3-none-any.whl` and `prepare-0.10.0.tar.gz`.

- [ ] **Step 2: Fresh-environment install**

```bash
python3.12 -m venv "$SCRATCHPAD/venv" 2>/dev/null || conda run -n metaquest python -m venv "$SCRATCHPAD/venv"
"$SCRATCHPAD/venv/bin/pip" install -q dist/prepare-0.10.0-py3-none-any.whl
"$SCRATCHPAD/venv/bin/prepare" --help | head -3
"$SCRATCHPAD/venv/bin/prepare" doctor; echo "exit=$?"
ls "$SCRATCHPAD/venv/bin" | grep -c metaquest
```

Expected: help names `prepare`; `doctor` runs (a non-zero exit for missing external tools is acceptable, an import error is not); the last count is 0.

- [ ] **Step 3: Merge**

```bash
git checkout main
git merge --ff-only rename/prepare
git branch -d rename/prepare
```

---

### Task 7: Outward-facing steps (main session, each confirmed with the user)

- [ ] **Step 1: Rename the GitHub repository**

```bash
gh repo rename prepare --repo FOI-Bioinformatics/metaquest --yes
git remote set-url origin git@github.com:FOI-Bioinformatics/prepare.git
git fetch origin && git status -sb | head -1
```

- [ ] **Step 2: Move the folder and copy the memory**

```bash
mv /Users/andreassjodin/Code/metaquest /Users/andreassjodin/Code/prepare
cp -R /Users/andreassjodin/.claude/projects/-Users-andreassjodin-Code-metaquest /Users/andreassjodin/.claude/projects/-Users-andreassjodin-Code-prepare
cd /Users/andreassjodin/Code/prepare && conda run -n metaquest pip install -e ".[dev]" -q && conda run -n metaquest prepare --help | head -1
```

- [ ] **Step 3: Push and tag**

```bash
cd /Users/andreassjodin/Code/prepare
git push origin main
git tag -a v0.10.0 -m "PREPARE 0.10.0: rename from MetaQuest"
git push origin v0.10.0
gh run list --repo FOI-Bioinformatics/prepare --limit 3
```

Expected: the CI and Release workflows start under the new names and finish green.
