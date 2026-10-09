# Rename MetaQuest to PREPARE

Date: 2026-10-09. Status: approved design.

## Purpose

The name MetaQuest clashes with a consumer product of the same name. The tool is renamed
PREPARE: the Python package, the command, the distribution, the GitHub repository, the
local folder, and every name the tool writes to disk. PREPARE is a plain name with no
expansion. The renamed release is 0.10.0.

## Decisions

- Package, command, PyPI distribution and conda environment are all `prepare`.
- Clean break on disk: no fallback to the old file names, environment variables or config
  path. The 0.10.0 CHANGELOG entry lists what a user of an existing project renames by hand.
- Version 0.10.0, not 1.0.0.
- Git history, past CHANGELOG entries, the historical specs and plans under
  `docs/superpowers/`, and the hidden old command names in `cli/commands/renamed.py` keep
  the old name.

## Name mapping

Applied case-aware across the tree outside the files listed as untouched.

| Old | New |
|---|---|
| `metaquest` (package, import path, command, PyPI, conda env, Makefile targets) | `prepare` |
| `MetaQuest`, `Metaquest` in prose and headings | `PREPARE` |
| `MetaQuest`, `Metaquest` in identifiers (`MetaQuestError`, `MetaquestParsePlugin`) | `Prepare` (`PrepareError`, `PrepareParsePlugin`) |
| `METAQUEST_<NAME>` settings variables, `METAQUEST_DATA`, `METAQUEST_PERF_SCALE` | `PREPARE_<NAME>` |
| `~/.config/metaquest/config.toml` | `~/.config/prepare/config.toml` |
| `metaquest_registry.json` | `prepare_registry.json` |
| `metaquest_store.json` (store catalogue and sidecar names), `metaquest_version` key | `prepare_store.json`, `prepare_version` |
| `.metaquest/runs/` | `.prepare/runs/` |
| `.metaquest-assembly.json` | `.prepare-assembly.json` |
| `fastq/.metaquest-tmp/` | `fastq/.prepare-tmp/` |
| `metaquest-lock-heartbeat` thread name | `prepare-lock-heartbeat` |
| sourmash entry points `metaquest_parse`, `_plot`, `_diversity`, `_taxonomy` | `prepare_parse`, `_plot`, `_diversity`, `_taxonomy` |
| `metaquest_plugins` entry-point group | `prepare_plugins` |
| logger names `metaquest.*` | `prepare.*` |
| `github.com/FOI-Bioinformatics/metaquest` | `github.com/FOI-Bioinformatics/prepare` |
| `~/Code/metaquest` | `~/Code/prepare` |
| GitHub Actions workflow display names `MetaQuest CI` and similar | `PREPARE CI` and similar |

Untouched: `.git/`, `CHANGELOG.md` entries before 0.10.0, `docs/superpowers/specs/*` and
`docs/superpowers/plans/*` dated before this spec, `renamed.py` old command names and their
tests, the `.claude/`, `.remember/` and `.superpowers/` folders. Memory notes under the
Claude memory folder are copied to the new path-keyed folder after the move, not edited.

## Approach

One mechanical commit first, then parallel fix-up by area, then a whole-branch review.

The mechanical commit, made on branch `rename/prepare` from main `de4a953`, does
`git mv metaquest prepare` and a scripted substitution in this order so that longer and
upper-case forms are replaced before shorter ones: `METAQUEST` to `PREPARE`, `MetaQuest`
and `Metaquest` to `PREPARE` in Markdown and to `Prepare` in Python, `metaquest` to
`prepare`. File names containing the old name are renamed with `git mv`. The script skips
the untouched paths. After that commit the tree imports and the suite runs.

The two stale worktrees under `.claude/worktrees/` (`reporting` and `storage-fixes`, both
clean and merged) are removed before any work, since their absolute paths break on the
folder move.

## Areas and owners

Four implementers on disjoint files after the mechanical commit, then one reviewer.

1. Package code (Opus): everything under `prepare/`. Checks each identifier the script
   produced, the plugin class names, the exception class, logger names, the thread name,
   settings prefix, config path, registry and store file names, run-log folder, assembly
   marker, scratch folder; and the matching unit tests for those names only
   (`tests/test_settings*.py`, `tests/test_registry*.py`, `tests/test_store*.py`,
   `tests/test_run_log*.py`, `tests/test_assembly*.py`).
2. Tests and tooling (Sonnet): every other file under `tests/`, `scripts/`, `Makefile`,
   `setup.cfg`, `pyproject.toml`, `environment.yml`, `requirements*.txt`, `local_test.sh`,
   `.github/workflows/*.yml`, `tests/test_module_sizes.py` KNOWN_EXCEPTIONS paths,
   `tests/helpers_processes.py`.
3. Documentation (Sonnet): `README.md`, `docs/*.md` (not `docs/superpowers/`),
   `CONTRIBUTING.md`, `CLAUDE.md`, `AGENTS.md` if tracked (it is not; it is edited in place
   but never added), the 0.10.0 entry in `CHANGELOG.md`, and the version bump in
   `pyproject.toml` and `prepare/__init__.py` (the only two lines of `pyproject.toml` this
   owner edits; owner 2 edits the rest).
4. Review (Opus): whole-branch review against this spec, then the residual search below.

## Residual search

After the fix-up, a case-insensitive search for `metaquest` over the tree, excluding
`.git/`, `.claude/`, `.remember/`, `.superpowers/`, `docs/superpowers/`, and the
pre-0.10.0 part of `CHANGELOG.md`, must return only `cli/commands/renamed.py` and the tests
of its hidden names, and the 0.10.0 CHANGELOG entry itself. `scripts/check_docs_commands.py`
and `scripts/check_no_print.sh` are updated so the quality gates look at `prepare/`.

## Verification

On the branch head, in order: `make check`, `make test`, `make pipeline`, `make build`;
then a fresh virtual environment installs the 0.10.0 wheel and runs `prepare --help` and
`prepare doctor`. `tests/test_optional_imports.py` must still pass with the renamed module
list. `tests/test_concurrency_processes.py` must pass since it spawns the real command and
reads `PREPARE_DATA`.

## Outward-facing steps

After the branch is merged to main and the gates are green, each step is confirmed with
the user before it runs:

1. `gh repo rename prepare` on `FOI-Bioinformatics/metaquest` (GitHub redirects the old URL).
2. `git remote set-url origin git@github.com:FOI-Bioinformatics/prepare.git`.
3. `mv ~/Code/metaquest ~/Code/prepare`, then copy the Claude memory folder
   `~/.claude/projects/-Users-andreassjodin-Code-metaquest/` to
   `~/.claude/projects/-Users-andreassjodin-Code-prepare/`.
4. `git push`, then tag and push `v0.10.0`; the release workflow builds the distribution.
   PyPI publishing stays gated on the `PYPI_PUBLISH` repository variable as before.

## Out of scope

PyPI publication, a bioconda recipe, any compatibility shim, the 0.9.1 left-overs listed in
the project memory, and any change to behaviour beyond the names.
