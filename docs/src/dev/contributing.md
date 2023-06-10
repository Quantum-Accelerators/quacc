# Contributing

## Steps to Contribute

To contribute to Quacc, we recommend doing the following:

1. [Fork the repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo) to your personal GitHub account.
2. [Clone this forked repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) to your local machine, e.g. via `git clone <repo url>.git`.
3. In the newly downloaded `quacc` base directory, run `pip install -e .[dev]` to install Quacc in editable mode and with the development dependencies.
4. [Commit your changes](https://github.com/git-guides/git-commit) and [push them](https://github.com/git-guides/git-push) to your personal forked repository.
5. Create a [pull request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to merge your changes into the main Quacc repository.

For reproducibility purposes, we strongly recommend installing Quacc in a fresh virtual environment.

## General Rules

Please abide by the following guidelines when contributing code to Quacc:

1. All changes should have associated unit tests.

2. All code should include type hints and have documentation for the inputs and outputs.

3. The input to most compute jobs should be an ASE `Atoms` object. The output of most compute tasks should be a schema from one of the module/functions within {obj}`quacc.schemas`.

4. Individual compute jobs should be written as functions, and workflows should be written as classes. Refer to {obj}`quacc.recipes.emt.core` and {obj}`quacc.recipes.emt.slabs` for examples.

5. Only define multi-step workflows if they go beyond simply stitching together existing functions.

6. When possible, you should use the "internal" geometry optimizers for a given code rather than the ASE optimizers.

7. Never use `@ct.lattice` in the code directly (unless it is wrapped by `@ct.electron` to make it a sublattice).

8. Ensure that the code remains flexible for the user.

9. `gzip` large test files to save space.

## Changelog

We keep a [`CHANGELOG.md`](https://github.com/quantum-accelerators/quacc/blob/main/CHANGELOG.md) file in the base directory of the `quacc` code. Before submitting your PR, be sure to update the `CHANGELOG.md` file under the "Unreleased" section with a brief description of your changes. The `CHANGELOG.md` file follows the [Keep a Changelog](https://keepachangelog.com) format.

## Style Guidelines

In general, please try to keep the code style consistent when possible. There are two main things to consider:

1. Functions should be lowercase and with underscores. Classes should be in CamelCase and constructed using a {obj}`@dataclass` decorator where possible.

2. All Python code should be formatted with [isort](https://github.com/PyCQA/isort) and then [black](https://github.com/psf/black), although this will be corrected automatically when merged.

## Unit Tests

All changes you make to Quacc should be accompanied by unit tests and should not break existing tests. To run the test suite, run `pytest` from the the `quacc/tests` directory.

If you are adding recipes based on a code that can be readily installed via `pip` or `conda` (e.g. tblite, xtb-python, DFTB+, Psi4), then you can run these codes directly in the test suite. Preferably, you should use a small molecule or solid and cheap method so the unit tests run quickly. If the recipes you're adding are proprietary or not available via `pip` or `conda` (e.g. Gaussian, GULP), then you will need to [monkeypatch](https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html) certain functions to change their behavior during testing. For instance, we do not want to run VASP directly during unit tests and have mocked the `atoms.get_potential_energy()` function to always return a dummy value of -1.0 during unit tests. Any mocked functions can be found in the `conftest.py` files of the testing directory.

## Tips

When developing new recipes, it is often helpful to start from an existing example. In general, we recommend referring to {obj}`quacc.recipes.emt` or {obj}`quacc.recipes.tblite`. Some advanced use cases can also be found in the VASP recipes.
