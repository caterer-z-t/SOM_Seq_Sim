# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- `som-fit --help` no longer crashes — removed an empty `metavar=""` on the
  `-t`/`-c`/`-o` arguments that triggered a CPython argparse formatting bug.
- `seq-sim` no longer crashes on every invocation — fixed an invalid
  `list(...)` call and a nonexistent `args.config` attribute reference in
  `Seq_Sim/seq_sim.py`.
- Corrected the README and paper Quick Start examples so the file names
  produced by `seq-sim` match the inputs expected by `som-fit`.
- Fixed a citation formatting typo and a help-text spacing typo in the SOM CLI.
- Corrected the LICENSE copyright holder list to match the paper authors.

### Added
- CI workflow to build and validate the distribution package (`build.yml`).
- CI workflow to publish releases to PyPI via GitHub Actions trusted
  publishing (`publish.yml`).

## [1.0.0] - 2024-12-02

### Added
- Command-line interfaces for both `Seq_Sim` (`seq-sim`) and `SOM` (`som-fit`).
- Hyperparameter tuning support in the SOM CLI (grid search over scaling
  method, grid dimensions, topology, and neighborhood function).
- Percent Variance Explained (PVE) and topographic error metrics for
  evaluating SOM fit quality.
- Component plane and categorical overlay plotting.
- `pytest` test suite with GitHub Actions CI and Codecov coverage reporting.
- Sphinx documentation hosted on Read the Docs, including tutorial notebooks
  (Iris, Titanic, and simulated sequencing data examples).
- Project logo, AI-usage disclosure, and JOSS submission materials
  (`paper.md`, packaging metadata).

### Changed
- Migrated example scripts from standalone `.py` files to Jupyter notebooks.
- Reorganized the repository into separate `SOM` and `Seq_Sim` packages.

## [0.0.0] - 2024-09-16

### Added
- Initial Self-Organizing Map implementation and pseudo-bulk sequencing data
  simulation.
- Initial GitHub Actions test workflow and `requirements.txt`.

[Unreleased]: https://github.com/caterer-z-t/SOM_Seq_Sim/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/caterer-z-t/SOM_Seq_Sim/compare/v0.0.0...v1.0.0
[0.0.0]: https://github.com/caterer-z-t/SOM_Seq_Sim/releases/tag/v0.0.0
