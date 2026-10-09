# Packaging and distribution

PREPARE ships as a source distribution and wheel on PyPI, built and attached to each GitHub
release by `.github/workflows/release.yml` (job `build`) whenever a `v*` tag is pushed. A second
job in that workflow, `pypi`, publishes the built distribution to PyPI; it stays off, skipped by
its own `if: vars.PYPI_PUBLISH == 'true'` condition, until a maintainer sets the repository
variable `PYPI_PUBLISH` to `true` (Settings > Secrets and variables > Actions > Variables) and the
"Before turning on PyPI publishing" checklist below is done. A bioconda recipe (a second, separate
distribution channel bioinformatics users often prefer) comes after PyPI publishing works, not
instead of it: bioconda's own review builds from a PyPI sdist, so there has to be one to point at
first. This repository does not carry a bioconda recipe folder; recipes for a bioconda-distributed
package live in the community's own `bioconda-recipes` repository, submitted as a pull request
there once PyPI publishing is live.

## Order of operations

1. Publish to PyPI (turn on `PYPI_PUBLISH`, tag a release, confirm `pip install prepare` works
   in a clean environment).
2. Submit a bioconda recipe (below) to `bioconda-recipes`, pointing its `source:` at that PyPI
   sdist.
3. Only after both are live, mention the conda install path in the README's install instructions
   alongside the existing pip one.

## Before turning on PyPI publishing

None of the following is checked by this change; they need a maintainer with network access and,
for the trusted-publisher step, a PyPI account with rights over the project (or the name).

- **PyPI trusted publisher**: `pypa/gh-action-pypi-publish` in the `pypi` job authenticates via
  PyPI's trusted-publisher (OIDC) mechanism, not a token, so it needs a trusted publisher
  registered on PyPI for this repository, workflow file (`release.yml`) and environment (`pypi`)
  *before* the first run. See <https://docs.pypi.org/trusted-publishers/>.
- **The name `prepare` on PyPI**: confirm at <https://pypi.org/project/prepare/> whether the
  distribution name is free. If it is taken, `pyproject.toml`'s `[project] name` becomes
  `prepare-bio` (the import name, `import prepare`, and the CLI command, `prepare`, both
  stay unchanged, only the *distribution* name `pip install prepare-bio` changes), and every
  reference to the distribution name in this document and the README's install instructions is
  updated to match.
- **`ncbi-datasets-cli` 16's flags**: `prepare/utils/tools.py`'s `TOOLS["datasets"]` pins a floor
  of 16 (`environment.yml`'s `ncbi-datasets-cli>=16`, checked by `tests/test_environment_pins.py`).
  Confirm that the `datasets download genome accession --inputfile ... --no-progressbar`
  invocation `prepare/core/constants.py`'s `ALLOWED_BIOINFORMATICS_TOOLS["datasets"]` allows
  still works unchanged on `ncbi-datasets-cli` 16.x before relying on it in a release announcement
  or a bioconda pin; that assumption is what set the floor at 16 rather than an older release.
- **bioconda `megahit` on `osx-arm64`**: confirm a `megahit` build exists for `osx-arm64` on the
  bioconda channel (`conda search -c bioconda megahit` on that platform, or
  <https://bioconda.github.io/recipes/megahit/README.html>) before telling Apple Silicon users
  conda install covers every tool PREPARE uses; the README already documents megahit 1.2.9's own
  Apple Silicon assembly-thread caveat (single-threaded by default on macOS), independent of
  whether the package itself is available.

## Bioconda recipe template

A starting point for the pull request to `bioconda-recipes`, once the checklist above is done.
Bioconda's own linter (`bioconda-utils lint`) and reviewers will ask for changes; this is a template, not
a recipe ready to submit as-is, and it is not stored anywhere else in this repository.

```yaml
{% set name = "prepare" %}
{% set version = "0.9.0" %}

package:
  name: "{{ name|lower }}"
  version: "{{ version }}"

source:
  url: "https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/prepare-{{ version }}.tar.gz"
  sha256: <sha256 of the sdist at that PyPI release; `pip download --no-binary :all: --no-deps
    prepare==0.9.0` then `sha256sum`, or the value PyPI's "Download files" page shows>

build:
  number: 0
  noarch: python
  script: "{{ PYTHON }} -m pip install . --no-deps --no-build-isolation -vv"
  entry_points:
    - prepare = prepare.cli.main:main

requirements:
  host:
    - python >=3.12
    - pip
    - setuptools >=61.0
    - wheel
  run:
    # Core runtime, from pyproject.toml [project] dependencies (the CLI must import and build
    # its parser with only these; see tests/test_optional_imports.py).
    - python >=3.12
    - pandas >=2.1
    - numpy >=1.26
    - matplotlib-base >=3.8
    - biopython >=1.82
    - lxml >=4.9.3
    - requests >=2.31
    # External tools, from prepare/utils/tools.py's TOOLS table (the same table
    # environment.yml is pinned against; see tests/test_environment_pins.py).
    - sra-tools >=3.0
    - minimap2 >=2.17
    - samtools >=1.10
    - megahit >=1.2.9
    - pigz >=2.4
    - ncbi-datasets-cli >=16

test:
  commands:
    - prepare --version
    - prepare doctor --json

about:
  home: "https://github.com/<org>/prepare"
  license: GPL-3.0-only
  license_file: LICENSE
  summary: "A toolkit for analyzing metagenomic datasets based on genome containment"
  doc_url: "https://github.com/<org>/prepare#readme"

extra:
  recipe-maintainers:
    - <github handle>
```

Notes on the template:

- The `run:` requirements list PREPARE's core (always-imported) runtime and the external tools
  from `TOOLS`; the optional extras (`analysis`, `interactive`, `maps`, `sourmash`, see
  `pyproject.toml`) are left out of this base recipe, the same way `pip install prepare` (no
  `[extra]`) leaves them out. A later, separate recipe or `run_constrained` block could offer
  them; that is a decision for whoever submits the bioconda pull request, not fixed here.
- `optional` tools in `TOOLS` (`prefetch`, `pigz`, `seqkit`, `megahit`) are still run
  dependencies here because bioconda recipes do not have PREPARE's own "falls back to a slower
  Python path" concept; the README's tool table documents which of them a user can skip and what
  happens then. `seqkit` is left out entirely, the same way it is commented out (not installed by
  default) in `environment.yml`.
- `matplotlib-base` (no Qt/GTK backend pulled in) is the usual bioconda choice for a
  headless CLI tool; confirm this still matches how PREPARE calls it
  (`prepare/visualization/`) before relying on it.
- The `sha256` and `home`/`recipe-maintainers` placeholders are filled in at submission time, once
  a real PyPI release and a bioconda-recipes contributor exist.
