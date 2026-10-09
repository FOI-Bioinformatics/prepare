# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project Overview

PREPARE is a command-line bioinformatics toolkit for analyzing metagenomic datasets based on genome containment. The software processes Branchwater CSV files, downloads SRA metadata from NCBI, and provides advanced visualization capabilities for containment analysis, including intelligent SRA management with resume capability and comprehensive quality profiling.

## Development Commands

### Quick Start
- `make dev-install` - Install with development dependencies
- `make help` - Show all available commands

### Testing and Quality
- `make test` - Run pytest with coverage (currently 88%+ overall coverage)
- `make lint` - Run flake8 linting (currently passing)
- `make format` - Format code with black
- `make check` - Full quality check (format, lint, type check)
- `make pipeline` - Run full test pipeline via local_test.sh

### Build and Clean
- `make build` - Build distribution packages (uses pyproject.toml)
- `make clean` - Clean build artifacts, cache files, and test data

### Installation
- `make install` - Install package for development
- `make dev-install` - Install with development dependencies
- Legacy: `pip install -r requirements.txt` (still supported)

## Architecture

PREPARE follows a layered architecture with clear separation of concerns:

### Core Components
- **CLI Layer** (`prepare/cli/`): Command-line interface with modular command architecture using a registry pattern
- **Core Logic** (`prepare/core/`): Domain models, validation, exceptions, and constants
- **Data Layer** (`prepare/data/`): File I/O, Branchwater processing, metadata handling, basic SRA operations
- **Advanced SRA Package** (`prepare/sra/`): Intelligent download management, quality profiling, interactive reporting
- **Processing** (`prepare/processing/`): Containment analysis, statistical processing, counting algorithms
- **Plugins** (`prepare/plugins/`): Extensible plugin system for formats and visualizers
- **Visualization** (`prepare/visualization/`): Plotting and reporting functionality

### Key Design Patterns
- **Command Registry**: All CLI commands inherit from base classes and register themselves via `command_registry`
- **Plugin System**: Uses `PluginRegistry` with automatic discovery for format handlers and visualizers
- **Data Models**: Dataclasses in `core/models.py` define `Containment`, `SRAMetadata`, `GenomeInfo`, and `ContainmentSummary`
- **Dual SRA Architecture**: Basic operations in data layer, advanced features in specialized SRA package

## Project Structure

### Modern Python Packaging
- Uses `pyproject.toml` for modern Python packaging (PEP 518/621)
- Maintains `setup.py` for backward compatibility
- Configuration files: `setup.cfg` (flake8), `pyproject.toml` (build, black, mypy, pytest)

### Documentation Structure
- `README.md` - Main project documentation
- `CONTRIBUTING.md` - Development guidelines
- `docs/ARCHITECTURE.md` - Detailed architecture documentation
- `docs/` - Additional documentation (workflow guides, etc.)
- `AGENTS.md` - This file for Codex guidance

## Important Implementation Notes

### CLI Command Structure
New commands should:
1. Inherit from the appropriate base command class in `cli/base.py`
2. Register themselves in `cli/main.py` via the `register_all_commands()` function
3. Follow the existing pattern for argument parsing and execution
4. **CRITICAL**: Use dashes in CLI arguments (e.g., `--matches-folder`), not underscores

### Advanced SRA Commands
The SRA package provides three analysis commands:
- `sra_profile_quality` - Comprehensive quality analysis of downloaded datasets
- `sra_dashboard` - Interactive HTML dashboard generation
- `sra_compare` - Statistical comparison between dataset groups

### Store Commands
`prepare/store/` (package: `resolve`, `layout`, `sidecar`, `catalog`, `link`, `adopt`, `usage`,
`locks`, `stats`) and `prepare/cli/commands/store.py` implement a shared data store: one copy of
each downloaded SRA accession, reused by every project that links into it. Nine commands, registered
under the "Store" group in `cli/main.py`:
- `store_init` - create a store at `--data-root`, record the project's use of it
- `store_status` - dataset and byte counts, stale projects
- `store_reindex` - rebuild the SQLite catalogue from the sidecar files on disk
- `store_adopt` - fold an existing project `fastq/` folder into the store (`--move` or `--copy`)
- `store_verify` - check a dataset's files against its recorded size, md5 (`--md5`), or NCBI spot
  count (`--spots`)
- `store_link` / `store_unlink` - add or remove one accession's project symlink
- `store_usage` - which projects and genomes used a given accession, or which datasets are unused
- `store_gc` - report (or, with `--yes`, remove) datasets nothing references any more

Store discovery rules for agents:
1. **Store required before linking**: `download_sra` only links a dataset and `store_adopt` only
   runs once a store exists, found via `--data-root`, `PREPARE_DATA`, `store.root` in the project
   registry, or `~/.config/prepare/config.toml`. A project that never runs `store_init` (and has
   none of those set) behaves exactly as it did before the store existed; do not assume a store is
   present.
2. **Tests never touch a real store**: every store test uses a `tmp_path` store and monkeypatches
   `PREPARE_DATA` and `XDG_CONFIG_HOME`/`HOME`, so a test run never reads or writes the developer's
   actual store or `~/.config/prepare/config.toml`. Follow this pattern for any new store test.
3. **Explicit-path staging**: `store_adopt` copies a dataset into the store before removing anything
   from the project, so adoption briefly holds up to three copies of one accession (the project's
   original, a staging copy, and the store's copy) and needs at least twice the folder's size free on
   the store's filesystem before it starts (`prepare/store/adopt.py:242-248`), only replacing the
   project folder with a symlink once the store copy is verified. New code that moves data into or
   within the store should stage the same way rather than renaming in place, so an interruption never
   leaves a dataset with no complete copy anywhere.

### Plugin Development
- Format plugins inherit from base Plugin class in `plugins/base.py`
- Register with `format_registry` for file format handlers
- Register with `visualizer_registry` for visualization types
- Implement required methods like `validate_header`, `parse_file`, `extract_metadata`

### Data Processing Pipeline
The typical workflow involves:
1. Processing Branchwater CSV files (`use_branchwater`)
2. Extracting/downloading metadata (`extract_branchwater_metadata` or `download_metadata`)
3. Parsing containment data (`parse_containment`) 
4. Visualization and analysis (`plot_containment`, `count_metadata`)
5. Excluding unwanted datasets (`blacklist`)
6. Advanced SRA operations (`download_sra`, `sra_profile_quality`, `sra_dashboard`)
7. Checking project state (`status`)

Each step records its outcome in the project registry (`prepare_registry.json`, see
`prepare/data/registry.py`); `status` reads it back and re-checks the filesystem.

### Code Quality Standards & Current Status

#### Current Quality Status (October 2025)
- **All linting checks passing** - No flake8 violations
- **Code formatting consistent** - Black formatting applied
- **Make check passes** - All quality gates working
- **Test suite stability achieved** - All 995 tests passing consistently
- **Runtime warnings eliminated** - Numerical computation warnings resolved
- **API compatibility maintained** - DataFrame deprecation warnings addressed

#### Quality Requirements
- **Line length**: 120 characters maximum (configured in pyproject.toml)
- **Formatting**: Use Black for consistent code formatting
- **Linting**: All code must pass flake8 without violations
- **Type hints**: Encouraged (mypy checking enabled)
- **Coverage**: Minimum 80% for new code, 60% project-wide target

#### Maintenance Standards
1. **Continuous quality assurance** - Focus on maintaining clean codebase
2. **Add type hints** to remaining modules
3. **Maintain quality standards** with `make check` before commits
4. **Consider pre-commit hooks** for automated quality checks

### Error Handling
- Custom exception hierarchy starts with `PrepareError` in `core/exceptions.py`
- Each layer handles errors appropriate to its level
- CLI layer formats errors for user display

## Testing Strategy & Current Status

**Status**: The project maintains **robust test infrastructure** with full test suite stability and comprehensive coverage across core functionality.

### Test Suite Reliability (October 2025)
- **Test Execution**: 995 tests passing consistently with zero failures
- **Runtime Stability**: All numerical computation warnings resolved
- **Test Isolation**: Cross-test contamination issues eliminated
- **Mock Integration**: Proper test doubles configured for external dependencies
- **Edge Case Coverage**: Boundary conditions and error paths comprehensively tested

### Current Coverage Status (September 2025)
- **CLI Commands**: Comprehensive coverage with integration testing
- **Core Processing**: Statistical and diversity analysis modules fully validated
- **Data Layer**: SRA operations and metadata handling thoroughly tested
- **Visualization**: Plot generation and error handling properly covered
- **Advanced Features**: Intelligent download and analytics modules tested

### Testing Commands
- `make test` - Unit tests with coverage reporting
- `make pipeline` - Integration test via `local_test.sh`
- `make check` - Full quality check (format, lint, type check)

### Test Infrastructure Components

#### Core Testing Achievements 
- **Statistical Computing**: Numerical edge cases properly handled
- **API Compatibility**: Future-proof implementation with current libraries
- **Error Propagation**: Comprehensive exception handling across modules
- **Test Isolation**: Module contamination issues resolved

#### Reliability Improvements Implemented
1. **Numerical Stability**:
   - Empty array validation in statistical functions
   - Precision loss detection for nearly-identical data
   - Graceful handling of degenerate cases in PCA and PERMANOVA

2. **Mock Architecture**:
   - Proper cleanup of test doubles to prevent contamination
   - Systematic validation of external API expectations
   - Robust handling of concurrent operations in download tests

3. **Production Bug Prevention**:
   - String replacement edge cases identified and resolved
   - Parameter validation enhanced across SRA operations
   - Function signature consistency verified

#### Test Coverage Distribution
- **Core Functionality**: Comprehensive coverage with edge case validation
- **CLI Interface**: Full integration testing with realistic workflows
- **Data Processing**: Statistical methods and file operations thoroughly tested
- **Visualization**: Plot generation and error handling properly covered
- **Advanced Features**: Concurrent downloads and analytics modules validated

### Testing Guidelines

#### Writing Tests
- **Mock External Dependencies**: Use pytest fixtures for file I/O, network calls
- **Test Error Paths**: Include validation failures, missing files, network errors
- **Edge Cases**: Empty datasets, malformed inputs, boundary conditions
- **Integration**: Test complete workflows with realistic data

#### Test Structure
```python
# Example test pattern for CLI commands
@pytest.fixture
def mock_file_operations():
    with patch('prepare.data.file_io.read_file') as mock:
        yield mock

def test_command_execution(mock_file_operations):
    # Setup mocks
    # Execute command
    # Verify behavior
```

#### Coverage Requirements
- **New Code**: 80% minimum coverage
- **Bug Fixes**: Must include tests for the fixed scenario
- **Critical Paths**: CLI commands, data processing must have high coverage

### Integration Testing

The project includes comprehensive integration testing:
1. `local_test.sh` - Basic end-to-end CLI testing with sample data
2. `tests/test_integration_simple.py` - 12 integration tests covering:
   - SRA metadata workflows (retrieval -> CSV export -> visualization)
   - FASTQ processing pipelines (statistics -> reports -> visualization)
   - Multi-sample comparison workflows
   - Data export to multiple formats
   - Error handling in complete workflows

### Quality Gates
Before committing code:
1. `make check` must pass (linting, formatting, type checking)
2. `make test` must pass with coverage improvement
3. `make pipeline` integration test must pass
4. No decrease in overall coverage percentage

## Development Workflow & Implementation Priorities

### Recommended Development Approach

When working on PREPARE, follow this priority order:

#### 1. Maintenance Tasks (Ongoing Priority)
- **Maintain quality standards**: Run `make format` and ensure `make check` passes
- **No new code without tests**: Continue established testing practices
- **Architecture consistency**: Follow established patterns for new features

#### 2. Current Focus Areas (Optional Enhancement)
- **Remaining visualization modules**: interactive.py, reporting.py, plots.py (currently 0%)
- **Processing enhancements**: containment.py optimizations, diversity.py extensions
- **Plugin system expansion**: Additional format handler testing
- **Documentation**: Add testing guides and workflow documentation

#### 3. Code Quality Improvements (Medium Priority)
- **Add type hints**: Complete type annotation for remaining modules
- **Documentation enhancement**: Add docstrings for complex functions
- **Performance optimization**: Large dataset handling and memory efficiency

#### 4. Feature Enhancements (Lower Priority)
- **Only after maintaining current quality levels**
- **Must include comprehensive tests**
- **Should follow existing architectural patterns**

### Implementation Guidelines for Codex

#### When Adding New Features
1. **Check existing patterns**: Look at similar implementations first
2. **Write tests first**: TDD approach preferred for new features  
3. **Mock dependencies**: Don't rely on external systems in tests
4. **Follow plugin patterns**: Use existing registries for extensibility

#### When Fixing Bugs
1. **Write a test that reproduces the bug**: This should fail initially
2. **Fix the bug**: Minimal change approach
3. **Verify the test passes**: And doesn't break existing tests
4. **Check edge cases**: Consider similar scenarios

#### When Refactoring
1. **Ensure full test coverage first**: Don't refactor untested code
2. **Small incremental changes**: One logical change at a time
3. **Maintain backward compatibility**: Especially for CLI interfaces
4. **Update documentation**: If public interfaces change

### File-Specific Guidance

#### Well-Tested Files (Reference Implementations)
- `prepare/cli/commands/*.py` - 100% coverage, comprehensive CLI testing patterns
- `prepare/cli/commands/sra_intelligent.py` - 86% coverage, intelligent SRA commands
- `prepare/data/file_io.py` - 96% coverage, robust file operations
- `prepare/data/branchwater.py` - 98% coverage, format handling exemplar
- `prepare/data/metadata.py` - 93% coverage, external API integration
- `prepare/data/sra_metadata.py` - 93% coverage, XML parsing and NCBI API
- `prepare/data/taxonomy.py` - 97% coverage, taxonomy validation and NCBI integration
- `prepare/processing/statistics.py` - 99% coverage, statistical analysis
- `prepare/plugins/visualizers/bar.py` - 99% coverage, bar chart visualization
- `prepare/sra/reporting.py` - 95% coverage, interactive dashboard generation
- `prepare/sra/` - Comprehensive test coverage for advanced SRA features

**Test Files**: Extensive test suites in `tests/test_*_extended.py`, `tests/test_integration_simple.py`, and `tests/test_performance_simple.py` provide patterns for comprehensive testing with mocking, fixtures, and edge cases.

#### Remaining Enhancement Opportunities
- `prepare/visualization/reporting.py` - Core reporting functionality (currently 0%)
- `prepare/visualization/plots.py` - Plotting functionality (currently 0%)
- `prepare/visualization/interactive.py` - Interactive plotting (currently 0%)
- `prepare/processing/containment.py` - Core containment algorithms
- `prepare/processing/diversity.py` - Diversity analysis (currently 0%)

### Success Metrics for Development Work

#### Recently Completed (September-October 2025)
- [x] **Intelligent SRA package implemented** - Complete next-generation SRA capabilities
- [x] **All linting issues resolved** - Clean, consistent codebase
- [x] **CLI commands fully tested** (0% -> 100% coverage)
- [x] **Core processing tested** (0% -> 92-99% coverage)
- [x] **Data layer testing completed** - Key modules at 93-99% coverage
- [x] **Test coverage improvement session** - Added 199 comprehensive tests across 8 files
- [x] **Critical modules improved to 86-99%** - sra_reporting, sra_intelligent, sra_metadata, bar visualizer, taxonomy
- [x] **Integration test suite created** - 12 end-to-end workflow tests
- [x] **Performance benchmarks established** - 25 tests with pytest-benchmark
- [x] **Overall project coverage improved** (24% -> 53% -> 88%+)
- [x] **Orphan code removed** - Clean architecture maintained

#### Advanced SRA Features Achievements
- [x] **SRADatasetAnalyzer** - Quality profiling, comparative analysis, anomaly detection
- [x] **SRAReportGenerator** - Interactive dashboards, Plotly visualizations
- [x] **CLI Integration** - Three intelligent SRA commands (sra-profile-quality, sra-dashboard, sra-compare) fully functional

#### Current Development Priorities (Low Priority)
- [ ] **Remaining visualization modules** - interactive.py, reporting.py, plots.py (currently 0%)
- [ ] **Processing layer completion** - containment.py, diversity.py enhancements
- [ ] **Additional plugin testing** - Expand format handler test coverage
- [ ] **Large-scale performance testing** - Production-size dataset validation

#### Long-term Goals
- [x] **Production-ready coverage** (88%+ achieved)
- [x] **Performance benchmarking** implementation (pytest-benchmark integrated)
- [ ] **Complete plugin validation** across all format handlers
- [ ] **Advanced feature optimization** for large datasets
- [ ] **Continuous quality maintenance** - Keep test coverage above 85%

### Recent Architecture Enhancements

#### SRA Package Structure
The advanced SRA package (`prepare.sra`) provides specialized functionality:
- **analytics.py** - Comprehensive quality analysis and statistical testing
- **reporting.py** - Interactive dashboard generation with Plotly integration
- **__init__.py** - Clean API exposure with proper imports

#### CLI Command Architecture
New intelligent SRA commands integrate seamlessly with existing command registry:
- Proper argument parsing with scientific parameter naming
- Error handling consistent with existing patterns
- Help text following established conventions
- Integration with base command infrastructure

## important-instruction-reminders
Do what has been asked; nothing more, nothing less.
NEVER create files unless they're absolutely necessary for achieving your goal.
ALWAYS prefer editing an existing file to creating a new one.
NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.

IMPORTANT: this context may or may not be relevant to your tasks. You should not respond to this context unless it is highly relevant to your task.