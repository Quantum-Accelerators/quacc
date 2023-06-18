# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

## [0.1.0]

### Added

- Added optional support for Covalent as the workflow manager.
- Improved documentation.
- Automatic formatting of files during PRs.
- This `CHANGELOG.md` file.
- The `BulkToSlabsFlow` in the EMT calculator for both Covalent and Jobflow.
- The `BulkToSlabsFlow` in the VASP calculator for Covalent.
- `LennardJones` demonstration recipes.
- A module named `quacc.util.db` to handle database interactions.
- A function `quacc.util.db.results_to_db` to store results in a database.
- A function `quacc.util.db.covalent_to_db` to store Covalent results in a database.
- A new `pop_analysis` kwarg to `quacc.schemas.cclib.summarize_run`.
- New functions `quacc.util.dict.clean_dicts` and `quacc.util.dicts.sort_dict`.
- Added support for Python 3.10.
- Added a name VASP setting `VASP_MIN_VERSION`.
- Added a new VASP swap for `EFERMI = "MIDGAP"` if `VASP_MIN_VERSION >= 6.4`.
- Monkeypatched in `.charge` and `.spin_multiplicity` properties for the `Atoms` object

### Changed

- Individual compute jobs are now functions instead of classes.
- All compute jobs are wrapped by `@ct.electron` instead of `@job`.
- Use of Jobflow is now optional.
- Use of a MongoDB is now optional.
- Changed the `extras` names in `setup.py`.
- Updated the dependencies.
- Increased default for `max_steps` in `run_ase_opt` from 100 to 500.
- Switch license to standard BSD-3.
- Simplified output schema for VASP slab workflows.
- Changed `removed_empties` default to `False` for the VASP schema.
- Use `setups_pbe54_MP` as the default pseudopotentials for the `BulkSet` preset in VASP.
- Modified thermo and vibrational frequency schemas.
- The `run_ase_opt` function returns an ASE `Optimizer` object instead of a `Trajectory` object.
- The `summarize_opt_run` function takes in an ASE `Optimizer` object instead of a `Trajectory` object.
- The `summarize_opt_run` function no longer needs the `atoms.calc.parameters` kwarg.
- Renamed `prerelax` kwarg in `qmof_relax_job` was changed to `run_prerelax` to reflect that it is a boolean.
- Renamed `volume_relax` kwarg in VASP recipes to `relax_cell` to reflect that it is a boolean.
- Renamed `molden` to `write_molden` in Gaussian recipes to reflect that it is a boolean.
- Renamed `quacc.recipes.tblite.core.thermo_job` to `quacc.recipes.tblite.core.freq_job`.
- Renamed `quacc.schemas.calc` to `quacc.schemas.ase`.
- Split `quacc.util.basics` and its associated functions into `quacc.util.dicts` and `quacc.util.files`.
- Moved `main` branch to an orphan and started new `main` branch.
- Moved ideal gas thermo calculation to `quacc.util.thermo`.
- Moved `check_logfile` to `quacc.util.files`.
- Moved `convert_auto_kpts` to `quacc.calculators.vasp`.
- By default, output schemas are sorted by key alphabetically.
- Used `|` operator for dictionary merging instead of `{**dict1, **dict2}`.
- Refactors VASP calculator.
- Added `relax_cell` option to EMT `relax_job`.

### Fixed

- If the user-specified `scratch_dir` does not exist on the filesystem when using `run_calc`, `run_ase_opt`, or `run_ase_vib`, it will be created.
- Fixed issue in the `TBLite` `ThermoJob` where magnetic moments were not being passed between the `Vibrations` and `IdealGasThermo` jobs.
- Fixed type hinting.

### Removed

- Removed the `xtb` recipes in favor of the `tblite` recipes following deprecation warning.
- Removed the `requirements.txt` file in favor of the `setup.py` file.
- Remove `quacc.util.dicts.merge_dicts`.
- Removed Python 3.8 support.

## [0.0.6]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.6 for more details.

## [0.0.5]

See https://github.com/quantum-accelerators/quacc/releases/tag/v0.0.5 for more details.

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
