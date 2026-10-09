# Contributing to PREPARE

Thank you for your interest in contributing to PREPARE! This document provides guidelines for contributing to the project.

## Development Setup

Requires Python 3.12 or newer.

1. **Clone the repository**:
   ```bash
   git clone https://github.com/FOI-Bioinformatics/PREPARE.git
   cd PREPARE
   ```

2. **Install in development mode**:
   ```bash
   pip install -e ".[dev]"
   ```
   Or, using conda: `make env`, `make env-dev` (installs into the `prepare` environment
   directly, no activation needed for this step), then `conda activate prepare` for day-to-day use.

3. **Run the test suite**:
   ```bash
   make test
   ```

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

## Development Workflow

1. **Code Quality**: Before submitting changes, ensure your code passes all quality checks:
   ```bash
   make check
   ```
   This runs:
   - `black` for code formatting
   - `flake8` for linting
   - `mypy` for type checking

2. **Auto-format code**:
   ```bash
   make format
   ```

3. **Run tests**:
   ```bash
   make test        # Run comprehensive test suite
   make pipeline    # Run full integration test
   ```

## Architecture Guidelines

### CLI Commands
- All commands should inherit from base classes in `cli/base.py`
- Register new commands in `cli/main.py`
- Follow existing argument naming conventions (use dashes, not underscores)

### Plugin Development
- Format plugins inherit from `Plugin` class in `plugins/base.py`
- Register with appropriate registries (`format_registry`, `visualizer_registry`)
- Implement required methods: `validate_header`, `parse_file`, `extract_metadata`

### Data Models
- Use dataclasses in `core/models.py` for data structures
- Follow the existing patterns for `Containment`, `SRAMetadata`, etc.

### Error Handling
- Use custom exceptions from `core/exceptions.py`
- Follow the hierarchy starting with `PrepareError`

## Testing

- Write tests for new functionality in the `tests/` directory
- Ensure both unit tests and integration tests pass
- The `local_test.sh` script provides end-to-end testing

### Network smoke test

`make test-network` downloads one 425-spot SRA run (SRR2517620) with fasterq-dump and
checks the files PREPARE writes. It needs `fasterq-dump` on the PATH and internet
access, so it is excluded from `make test`. Run it before releasing any change to
the `prepare/data/sra/` package or `prepare/utils/security.py`.

## Submitting Changes

1. **Create a feature branch**: `git checkout -b feature/your-feature-name`
2. **Make your changes** following the guidelines above
3. **Run quality checks**: `make check`
4. **Run tests**: `make test`
5. **Commit your changes**: Use clear, descriptive commit messages
6. **Push and create a pull request**

## Code Style

- Follow PEP 8 Python style guidelines
- Use `black` for automatic code formatting (line length: 88)
- Use type hints where appropriate
- Write clear docstrings for functions and classes

## Testing Guidelines

- **Comprehensive Coverage**: Include unit tests for all new functionality
- **Edge Case Handling**: Test boundary conditions and degenerate cases
- **Numerical Stability**: Verify proper handling of empty arrays, near-zero variances, and precision loss
- **Mock External Dependencies**: Use proper test doubles for NCBI APIs, file operations, and system calls
- **Test Isolation**: Ensure tests do not contaminate each other's state
- **Error Path Testing**: Include tests for validation failures and error conditions

## Pull Request Guidelines

- Provide a clear description of the changes
- Include tests for new functionality with proper edge case handling
- Ensure all tests pass without warnings
- Maintain numerical stability in statistical computations
- Update documentation if necessary
- Reference any related issues

## Questions?

If you have questions about contributing, please open an issue or start a discussion in the repository.