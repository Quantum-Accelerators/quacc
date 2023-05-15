# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

### Added

### Changed

### Fixed

## [0.1.0]

See https://github.com/arosen93/quacc/releases/tag/v0.1.0 for more details.

### Added
- Added support for Covalent as the workflow manager.
- Lots of new documentation.
- Type checking during PRs.
- Automatic formatting of files during PRs.
- This `CHANGELOG.md` file.

### Changed
- Individual compute jobs are now functions instead of classes.
- There is no longer an explicit dependency on Jobflow, relying on the user to specify the `@job` and `Flow()` wrappers.
- Use of a MongoDB is now optional to make setup easier.
- Split `quacc.util.basics` into `quacc.util.dicts` and `quacc.util.files`.
- Changed the `extras` names in `setup.py`.

### Removed
- Removed the `xtb` recipes in favor of the `tblite` recipes.
- Removed the `requirements.txt` file in favor of the `setup.py` file.

### Fixed

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