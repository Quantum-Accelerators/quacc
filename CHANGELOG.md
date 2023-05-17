# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

### Added

- Added optional support for Covalent as the workflow manager.
- Improved documentation.
- Automatic formatting of files during PRs.
- This `CHANGELOG.md` file.
- The `BulkToSlabsFlow` in the EMT calculator for both Covalent and Jobflow.
- The `BulkToSlabsFlow` in the VASP calculator for Covalent.
- `LennardJones` demonstration recipes.

### Changed

- Individual compute jobs are now functions instead of classes.
- All compute jobs are wrapped by `@ct.electron` instead of `@job`.
- Use of Jobflow is now optional.
- Use of a MongoDB is now optional.
- Split `quacc.util.basics` and its associated functions into `quacc.util.dicts` and `quacc.util.files`.
- Changed the `extras` names in `setup.py`.
- Updated the dependencies.
- Increased default for `max_steps` in `run_ase_opt` from 100 to 500.
- Renamed `prerelax` kwarg in `qmof_relax_job` was changed to `run_prerelax` to reflect that it is a boolean.
- Renamed `volume_relax` kwarg in VASP recipes to `relax_volume` to reflect that it is a boolean.
- Renamed `molden` to `write_molden` in Gaussian recipes to reflect that it is a boolean.
- Switch license to standard BSD-3.
- Simplified output schema for VASP slab workflows.
- Fixed type hinting.

### Fixed

- If the user-specified `scratch_dir` does not exist on the filesystem when using `run_calc`, `run_ase_opt`, or `run_ase_vib`, it will be created.

### Removed

- Removed the `xtb` recipes in favor of the `tblite` recipes following deprecation warning.
- Removed the `requirements.txt` file in favor of the `setup.py` file.

## [0.0.6]

See https://github.com/arosen93/quacc/releases/tag/v0.0.6 for more details.

## [0.0.5]

See https://github.com/arosen93/quacc/releases/tag/v0.0.5 for more details.

## [0.0.5]

See https://github.com/arosen93/quacc/releases/tag/v0.0.5 for more details.

## [0.0.4]

See https://github.com/arosen93/quacc/releases/tag/v0.0.4 for more details.

## [0.0.3]

See https://github.com/arosen93/quacc/releases/tag/v0.0.3 for more details.

## [0.0.2]

See https://github.com/arosen93/quacc/releases/tag/v0.0.2 for more details.

## [0.0.1]

See https://github.com/arosen93/quacc/releases/tag/v0.0.1 for more details.
