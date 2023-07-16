# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

### Added

- Added a custom Q-Chem calculator based on Pymatgen and Custodian.
- Added Q-chem core recipes.
- Added the ability to pass option `run_kwargs` to the `.run()` method of the ASE optimizers.
- Added a `PRIMARY_STORE` quacc setting that, when specified, will automatically store all results in the specified Maggma store.

### Changed

- The trajectory is now stored in `quacc.schemas.cclib.summarize_run`.

### Fixed

- Psi4 recipes now use `uks` for unrestricted calculations instead of `uhf`.

## [0.1.2]

### Added

- Added support for Prefect.
- Added `quacc.util.dask` module.

### Changed

- Made slab-related kwargs more consistent, such as by changing `slab_relax_electron`/`slab_relax_app` to just `slab_relax` regardless of workflow manager.

## [0.1.1]

### Added

- Added this `CHANGELOG.md` file.
- Added the `quacc.recipes.emt.parsl` module.
- Added a CLI for `quacc` along with a `quacc config` option to configure Covalent appropriately upon install.
- Added generic type hints for schemas.
- Added a `CREATE_UNIQUE_WORKDIR` global setting to have quacc automatically make a unique working directory for each calculation.
- Added `CHECK_CONVERGENCE` to global settings.

### Changed

- The `quacc` directory is now found in a `src` basefolder.
- All recipes now support the `Atoms` object being passed in as an `AtomsSchema`.
- The `slab_relax_job` kwarg in `recipes.emt.slabs` and `recipes.vasp.slabs` workflows can no longer be set to `None`, as there are few situations where this would be desired.
- Class-based recipes have been converted to functions since they don't save state or have inheritance.

### Docs

- Switched the docs from Furo to Material for MkDocs.
- Expanded upon Parsl documentation.
- Modified tutorials and added example configs for Slurm.

### Fixed

- Fixed a bug where parallel, multithreaded Python jobs could cause file I/O issues.
- Temporary directories are cleaned up after the run is completed.

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
