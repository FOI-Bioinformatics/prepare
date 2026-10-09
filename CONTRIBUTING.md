# Contributing to PREPARE

This document describes how to set up a development environment, check a change and submit it.

## Development setup

Requires Python 3.12 or newer.

```bash
git clone https://github.com/FOI-Bioinformatics/prepare.git
cd prepare
pip install -e ".[dev]"
make test
```

With conda, `make env` creates the `prepare` environment from `environment.yml` (Python 3.12 and the
external tools), and `make env-dev` installs the package with its development dependencies into it.
`make env-dev` runs `conda run -n prepare pip install -e ".[dev]"`, so it installs into that
environment whatever environment is active. Run `conda activate prepare` afterwards to use the
`prepare` console script day to day. Do not keep PREPARE installed in another interpreter's
site-packages at the same time; run `pip uninstall prepare` there, so the console script resolves to
the copy in the `prepare` environment.

### Several worktrees and one environment

An editable install (`pip install -e .`) points the environment at the one checkout it was run
from. With several git worktrees sharing the `prepare` environment, `python -m pytest` or
`python -c "import prepare"` run from a worktree's root still imports that worktree's code,
because the current directory comes first on `sys.path`. A script run from anywhere else, or the
`prepare` console script, imports whichever checkout the editable install points at, which may
be a different branch. This is why `scripts/smoke_chain.sh` sets `PYTHONPATH` to its own
repository root before it changes into the scratch project directory. Do the same
(`PYTHONPATH=/path/to/worktree`) when running a helper script from outside a worktree, or
re-run `pip install -e .` in the worktree you are working in.

## Checking a change

```bash
make format      # black, line length 120
make check       # every quality gate
make test        # the test suite with a coverage report
make pipeline    # end-to-end run of local_test.sh
make help        # every make target
```

`make check` runs, in order: black, flake8, mypy, the no-print check, the plotly CDN check, the radon
complexity ceiling, the module size and maintainability ceilings, the ASCII check, the atomic-write
check, and the check that README.md names every visible command and no unknown one. The `check`
target in the Makefile is the authoritative list, and the "Quality gates" section of
[CLAUDE.md](CLAUDE.md) explains each gate. `pip-audit --strict` runs in CI, not in `make check`.

### Network smoke test

`make test-network` downloads one 425-spot SRA run (SRR2517620) with fasterq-dump and
checks the files PREPARE writes. It needs `fasterq-dump` on the PATH and internet
access, so it is excluded from `make test`. Run it before releasing any change to
the `prepare/data/sra/` package or `prepare/utils/security.py`. The same chain runs nightly in CI
(`.github/workflows/smoke.yml`): it downloads the run, validates and profiles it, extracts reads against
the test genome and writes the results table. It is the only CI job that uses the network.

## Code guidelines

- **Commands** inherit from the base classes in `prepare/cli/base.py`, register through
  `register_all_commands()` in `prepare/cli/main.py`, and use dashes in flag names
  (`--matches-folder`). [CLAUDE.md](CLAUDE.md) lists the further rules a new command follows
  (output through `emit`, run log, required tools, settings, errors).
- **Plugins** subclass `Plugin` in `prepare/plugins/base.py`, which defines only a name, a
  description and a version, and register with `format_registry` or `visualizer_registry` from the
  same module. The methods a format handler provides follow the Branchwater handler in
  `prepare/plugins/formats/branchwater.py`. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).
- **Data models** are the dataclasses in `prepare/core/models.py`.
- **Errors** use the exceptions in `prepare/core/exceptions.py`, whose hierarchy starts with
  `PrepareError`; each class carries its exit code.
- **Style** is PEP 8 as enforced by black (line length 120) and flake8, with type hints and
  docstrings on public functions and classes.

## Tests

- Add tests for new behaviour in `tests/`, including error paths and boundary cases.
- A bug fix comes with a test that fails before the fix.
- Mock external services (NCBI, Branchwater, external tools); only `make test-network` uses the
  network.
- Keep tests independent of each other's state, and of the developer's own store and
  configuration file.

## Submitting changes

1. Create a branch: `git checkout -b feature/your-feature-name`.
2. Make the change, with tests and documentation where needed.
3. Run `make check` and `make test`.
4. Commit with a clear message, push, and open a pull request that describes the change and refers to
   any related issue.

Questions can be raised as an issue in the repository.
