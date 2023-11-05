# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0]

### Added

- Added `quacc.recipes.common`, including `quacc.recipes.common.core` and `quacc.recipes.common.phonons`
- Added phonopy schema in `quacc.schemas.phonopy`
- Added EMT phonon workflow

### Changed

- The method for specifying swaps to the default calculator arguments in recipes has been overhauled. In place of the previous `calc_swaps: dict` approach, additional parameters can be specified via the `**kwargs`. In practice, this means in recipes you can do `scf = "diis", max_steps=100` instead of `calc_swaps={"scf":"diis", "max_steps":100}`, for instance.
- Moved `quacc.runners.calc` to `quacc.runners.ase`
- Moved `quacc.runners.prep` to `quacc.schemas.prep`

## [0.3.13]

### Fixed

- Fixed support for the most recent version of ASE as of 11/2/2023

## [0.3.12]

### Added

- Added more type hints for input keyword arguments
- Added the `quacc info` command-line argument for error reporting

### Changed

- Consolidated VASP INCAR co-pilot settings
- Internal refactoring of the Q-Chem and VASP calculators

### Fixed

- All `Path`-based settings have `.expanduser()` applied

## [0.3.11]

### Changed

- Improved error handling

## [0.3.10]

### Changed

- DFTB+ `MaxSccIterations` set to 200 by default
- SCC errors in DFTB+ relaxations will be reported with higher priority than a relaxation error

## [0.3.9]

### Changed

- Changed prefix of scratch directories to have the datetime in it for easier debugging

### Fixed

- DFTB+ recipes now won't crash when `kpts=None` for structures with PBCs

## [0.3.8]

### Changed

- Changed minimum version of Parsl

## [0.3.7]

### Added

- Added type-hinting for schemas

### Changed

- Recipes now take `Atoms` instead of `Atoms | dict` as input
- Recipes no longer require the use of `fetch_atoms`
- Standardized schema fields

### Removed

- Removed the `fetch_atoms` function

## [0.3.6]

### Fixed

- Fixed VASP recipes to always respect `calc_swaps`.
- Fixed setting multiple cores in Custodian for Q-Chem runs.

## [0.3.5]

### Added

- Added a new `VibThermoSchema` that merges the `VibSchema` and `ThermoSchema` as a convenience

### Changed

- The output of recies no longer automatically removes blank lists/dicts/sets to maintain a consistent schema.
- The tblite `freq_job` and lj `freq_job` now use `VibThermoSchema` as the output
- Changed default `SCRATCH_DIR` from `Path.cwd() / ".scratch"` to `Path("~/.scratch")`

## [0.3.4]

### Added

- Added a new documentation section for worked examples on HPC.

### Changed

- Updated the minimum required versions for Pymatgen and emmet-core to reflect Pydantic 2 update.
- Added `symprec=1e-8` to the default parameters of VASP relaxation jobs.
- VASP bulk recipes now use "BulkSet" as the default preset instead of `None`.
- VASP slab recipes now use "SlabSet" as the default preset instead of `None`.

## [0.3.3]

### Fixed

- Fixed pip version conflict when installing Covalent dependencies.

## [0.3.2]

### Changed

- Migration to Pydantic 2 is complete, along with the corresponding dependencies.

## [0.3.1]

### Changed

- Pydantic field type validation has been improved.
- Moved `quacc.utils.atoms`/`.slabs`/`.defects` into `quacc.atoms.core`/`.slabs`/`.defects`
- Moved `quacc.custodian` to `quacc.calculators.custodian`
- Moved `quacc.utils.thermo.ideal_gas` to `quacc.builders.thermo.build_ideal_gas`
- Renamed `quacc.schemas.ase.summarize_thermo` to `quacc.schemas.ase.summarize_ideal_gas_thermo`

## [0.3.0]

### Changed

- Changed default `SCRATCH_DIR` from `/tmp` to `Path.cwd() / .scratch`
- Refactored decorator handling to be more concise
- Renamed `[optimizers]` extras to `[sella]`
- Moved `quacc.utils.wflows` into `quacc.wflow.decorators` and `quacc.wflow.prefect`
- Moved `quacc.utils.db` into `quacc.wflow.db`
- Moved `quacc.utils.calc` to `quacc.runners.ase`
- Moved `quacc.presets` to `quacc.calculators.presets`

### Fixed

- Fixed parsing of `RESULTS_DIR` and `SCRATCH_DIR` paths when they contain `~`
- Symlinks are no longer created if `RESULTS_DIR == SCRATCH_DIR` since this is redundant
- If `RESULTS_DIR` or `SCRATCH_DIR` are not present, they will be automatically made

## [0.2.8]

### Changed

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
- Added `quacc.wflow` module.

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
