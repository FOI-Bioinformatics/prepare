.PHONY: env env-dev test test-perf test-network lint format check doctor build clean pipeline install help dev-install

help:
	@echo "Available commands:"
	@echo "  env         - Create or update the 'prepare' conda environment (Python 3.12)"
	@echo "  env-dev     - Install the package with dev extras into the 'prepare' conda environment"
	@echo "  test        - Run tests with coverage"
	@echo "  test-perf   - Run only the performance bounds (scale with PREPARE_PERF_SCALE)"
	@echo "  lint        - Run flake8 linting"
	@echo "  format      - Format code with black"
	@echo "  check       - Run all quality checks (format, lint, type check, module size, ASCII, atomic writes, docs)"
	@echo "  doctor      - Check external tools, configuration, store, disk space and resources"
	@echo "  build       - Build distribution packages"
	@echo "  clean       - Clean build artifacts and cache"
	@echo "  pipeline    - Run full integration test pipeline"
	@echo "  install     - Install package for development"
	@echo "  dev-install - Install with development dependencies"

env:
	conda env create -f environment.yml || conda env update -f environment.yml --prune

env-dev:
	conda run -n prepare pip install -e ".[dev]"

dev-install:
	python -m pip install -e ".[dev]"

install:
	python -m pip install -e .

test:
	python -m pytest tests/ --cov=prepare

test-perf:
	python -m pytest -m perf -q

test-network:
	python -m pytest tests/test_network_smoke.py -m network -x -v

lint:
	python -m flake8 prepare tests

format:
	python -m black prepare tests

check:
	@echo "Running under $$(python --version)..."
	@echo "Running format check..."
	python -m black --check --diff prepare tests
	@echo "Running linting..."
	python -m flake8 prepare tests
	@echo "Running type check..."
	python -m mypy prepare
	@echo "Checking that nothing prints outside prepare/cli/base.py (commands use self.emit)..."
	bash scripts/check_no_print.sh
	@echo "Guarding against the frozen plotly-latest CDN alias..."
	@if grep -rn "cdn.plot.ly/plotly-latest" prepare --include='*.py'; then \
		echo "ERROR: use prepare.utils.html.plotly_cdn_script() instead of the plotly-latest alias"; \
		exit 1; \
	fi
	@echo "Checking cyclomatic complexity ceiling (fail on rank D or worse)..."
	@output=$$(python -m radon cc prepare -n D -s); \
	if [ -n "$$output" ]; then \
		echo "$$output"; \
		echo "ERROR: functions at complexity rank D or worse; refactor before merging"; \
		exit 1; \
	fi
	@echo "Checking module size (800 lines) and maintainability (MI >= 20) ceilings..."
	python tests/test_module_sizes.py --check
	@echo "Checking for non-ASCII bytes outside the documented exemptions..."
	bash scripts/check_ascii.sh
	@echo "Checking that output files are written atomically (temp name, then rename)..."
	bash scripts/check_atomic_writes.sh
	@echo "Checking every command is documented and README names no unknown command..."
	python scripts/check_docs_commands.py
	@echo "All quality checks passed!"

doctor:
	python -m prepare.cli.main doctor

build:
	python -m build

clean:
	rm -rf build/ dist/ *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -name ".coverage" -delete

pipeline:
	bash local_test.sh
