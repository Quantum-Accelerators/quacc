# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

### Added

- Added this `CHANGELOG.md` file.
- Added the `quacc.recipes.emt.parsl` module.
- Added a CLI for `quacc` along with a `quacc config` option to configure Covalent appropriately upon install.
- Added a `quacc.util.wflows.make_dask_cluster` function to easily generate a Dask cluster.
- Added generic type hints for schemas.
- Added `CHECK_CONVERGENCE` to global settings.
- Added a custom Q-Chem calculator based on Pymatgen and Custodian.
- Added Q-chem core recipes.

### Changed

- All recipes now support the `Atoms` object being passed in as an `AtomsSchema`.
- The `slab_relax_job` kwarg in `recipes.emt.slabs` and `recipes.vasp.slabs` workflows can no longer be set to `None`, as there are few situations where this would be desired.
- Class-based recipes have been converted to functions since they don't save state or have inheritance.

### Docs

- Expanded upon Parsl documentation.
- Modified tutorials and added example configs for Slurm.

### Fixed

- Fixed a bug where parallel, multithreaded Python jobs could cause file I/O issues.

### Removed

- Removed `quacc.recipes.vasp.jobflow` module to prioritize Covalent and Parsl.

## [0.1.0]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.1.0 for more details.

## [0.0.6]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.6 for more details.

## [0.0.5]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.5 for more details.

## [0.0.4]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.4 for more details.

## [0.0.3]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.3 for more details.

## [0.0.2]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.2 for more details.

## [0.0.1]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.1 for more details.
