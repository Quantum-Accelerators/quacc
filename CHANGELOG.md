# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.8]

###

- Decorators are now applied at import time rather than at function-call time.

### Fixed

- Fixed compatability with FireWorks.
- Fixed I/O issue on Windows and NFS filesystems where Sella optimizations would crash.

## [0.2.7]

### Added

- Added a `freq_job` recipe in Q-Chem (@samblau)
- Added support for TRICs in Sella calculations (@samblau)

### Changed

- Refactored all recipes to prevent code duplication
- Set `LREAL = False` in VASP `BulkSet.yaml`

## [0.2.6]

### Added

- Add support for `Path` objects everywhere possible.

### Changed

- Charge and spin multiplicity are now required arguments in molecular DFT calculators.
- Slab recipes now use `make_slabs_from_bulk` instead of `make_max_slabs_from_bulk`
- Use the `logging` module when warnings do not need to be immediately addressed.
- Functions are no longer used as kwargs in recipes to help with (de)serialization in certain workflow engines.

### Removed

- Removed `make_max_slabs_from_bulk` function.

## [0.2.5]

### Added

- Added support for the Redun workflow engine.
- Added a simple CLI for (un)setting YAML config parameters.

### Changed

- Internal refactoring of the `Vasp` calculator.

### Removed

- Removed the demonstration Jobflow EMT slab recipe.

## [0.2.4]

### Changed

- Updated VASP co-pilot rules and `BulkSet.yaml` in accordance with recent MP benchmarking (#849).
- Moved `fetch_atoms` to `quacc.schemas`.

### Fixed

- VASP: Fixed a scenario where SIGMA was set to 0.05 if the user provided SIGMA > 0.05. This was supposed to happen when ISMEAR = 0 but was previously done when ISMEAR = -5, which is not influenced by SIGMA.
- Fixed `mp_relax_flow` to use a `#!Python @flow` decorator instead of `#!Python @job`
- Fixed the Q-Chem Quasi-IRC job to run as one `#!Python @job` instead of two.

## [0.2.3]

### Added

- Added new defect EMT workflows (@rwexler)
- Added new defect utilities (@rwexler)
- Added a custom Q-Chem calculator based on Pymatgen and Custodian. (@samblau)
- Added Q-chem core recipes. (@samblau)
- Added the ability to pass option `run_kwargs` to the `.run()` method of the ASE optimizers.

### Changed

- Removed the need to call `ct.dispatch` when launching Covalent workflows.
- Renamed `quacc.util` to `quacc.utils`
- Refactored and standardized the NewtonNet recipes to match the other recipes

### Fixed

- The `WORKFLOW_ENGINE` setting can now be set dynamically, as expected.

## [0.2.2]

### Changed

- Changed package handling in `pyproject.toml` to better handle pip installs.

## [0.2.1]

### Added

- Added NewtonNet recipes and better sella support (@kumaranu)

### Changed

- Standardized `relax_cell` to be `False` by default for all recipes.
- Update the default `SIGMA` and `ISMEAR` in `MPScanSet` to match the new values to be used in Atomate2.

### Fixed

- Fixed a bug on Windows where the trajectory file would not be closed after a relaxation, causing permission errors.

### Removed

- Removed atomate2 as a dependency.
- Removed experimental Prefect support.

## [0.2.0]

### Added

- Added a quacc global setting, `WORKFLOW_ENGINE`, to set the workflow manager.
- Added support for the new `covalent-hpc-plugin` in `pyproject.toml`.

### Changed

- Quacc recipes are now decorated with generic `#!Python @job`, `#!Python @subflow`, or `#!Python @flow` decorators to be workflow engine-agnostic.
- Covalent is now an optional dependency to allow the base `quacc` package to be dependency-light.
- Default VASP pseudopotentials changed to `setups_pbe54.yaml`.
- Changed `quacc.recipes.emt.jobflow` and `quacc.recipes.emt.prefect` to `quacc.recipes.emt._jobflow` and `quacc.recipes.emt._prefect`

### Fixed

- Set `LMAXMIX` based on the `Z` value, not the s/p/d/f-block type.
- Do not set automatically `LMAXTAU` to 8 for f-containing elements in `Vasp` calculator.

### Removed

- Removed `quacc.recipes.emt.parsl` since it is now obsolete.
- Removed `quacc config` command-line option.

## [0.1.3]

### Added

- Added a `PRIMARY_STORE` quacc setting that, when specified, will automatically store all results in the specified Maggma store.

### Changed

- The trajectory is now stored in `quacc.schemas.cclib.summarize_run`.
- By default, results are now stored in a fixed `RESULTS_DIR`, which defaults to the current working directory.

### Fixed

- Fix path-related issues when running in local multi-threaded mode.
- Psi4 recipes now use `uks` for unrestricted calculations instead of `uhf`.

## [0.1.2]

### Added

- Added support for Prefect.
- Added `quacc.utils.wflows` module.

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
