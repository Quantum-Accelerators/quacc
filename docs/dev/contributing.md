# Contributing

We welcome all contributions, big or small (provided they come with unit tests!). Here are some general guidelines and tips. Note that it is never too early to open a PR if you're looking for input or feedback!

## General Tips

When developing new recipes, it is often helpful to start from an existing example. In general, we recommend referring to `quacc.recipes.emt` or `quacc.recipes.tblite`. Some advanced use cases can also be found in the VASP recipes.

## Steps to Contribute

To contribute to quacc, we recommend doing the following:

- [Fork the repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo) to your personal GitHub account.

- [Clone this forked repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) to your local machine, e.g. via `git clone <repo url>.git`.

- In the newly downloaded `quacc` base directory, run `pip install -e .[dev]` to install quacc in editable mode and with the development dependencies.

- [Commit your changes](https://github.com/git-guides/git-commit) and [push them](https://github.com/git-guides/git-push) to your personal forked repository _in a new branch_.

- Create a [pull request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to merge your changes into the main quacc repository.

For reproducibility purposes, we strongly recommend installing quacc in a fresh virtual environment.

## Guidelines

Please abide by the following guidelines when contributing code to quacc:

- All changes should have associated unit tests that pass. Code coverage should be maintained.

- All code should include type hints and have internally consistent documentation for the inputs and outputs.

- Only define multi-step workflows if they go beyond simply stitching together existing functions or if they are widely used in other recipes. Otherwise, just define the individual functions.

- Ensure that the code remains flexible for the user whenever possible.

- Where appropriate, you should use the "internal" geometry optimizers for a given code rather than the ASE optimizers.

- Use absolute paths where possible and avoid using `os.chdir`. This helps ensure that quacc can be run in a multi-threading mode.

- `gzip` large test files to save space.

- Update the `CHANGELOG.md` file.

## Changelog

We keep a `CHANGELOG.md` file in the base directory of the quacc code. Before submitting your PR, be sure to update the `CHANGELOG.md` file under the "Unreleased" section with a brief description of your changes. The `CHANGELOG.md` file follows the [Keep a Changelog](https://keepachangelog.com) format.

## Style Guidelines

In general, please try to keep the code style consistent when possible. There are two main things to consider:

1. All Python code should be formatted with [isort](https://github.com/PyCQA/isort) and then [black](https://github.com/psf/black), although this will be corrected automatically when merged.

2. Use NumPy-style Docstrings.

3. Address any relevant style issues raised by the "Deepsource: Python" and "Sourcery" GitHub actions, if present.

4. Functions should be lowercase and with underscores. Classes should be in PascalCase and constructed using a `#!Python @dataclass` decorator where possible.

## Unit Tests

All changes you make to quacc should be accompanied by unit tests and should not break existing tests. To run the full test suite, run `pytest .` from the the `quacc/tests` directory. Each PR will report the coverage once your tests pass, but if you'd like to generate a coverage report locally, you can use [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/), such as by doing `pytest --cov=quacc .` in the `tests` directory.

If you are adding recipes based on a code that can be readily installed via `pip` or `conda` (e.g. tblite, DFTB+, Psi4), then you can run these codes directly in the test suite. Preferably, you should use a small molecule or solid and cheap method so the unit tests run quickly. If the recipes you're adding are proprietary or not available via `pip` or `conda` (e.g. Gaussian, GULP), then you will need to [monkeypatch](https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html) certain functions to change their behavior during testing. For instance, we do not want to run VASP directly during unit tests and have mocked the `atoms.get_potential_energy()` function to always return a dummy value of -1.0 during unit tests. Any mocked functions can be found in the `conftest.py` files of the testing directory.
