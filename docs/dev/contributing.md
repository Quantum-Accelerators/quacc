# Contributing Guide

We welcome all contributions, big or small (provided they come with unit tests!). Here are some general guidelines and tips. Note that it is never too early to open a PR if you're looking for input or feedback!

!!! Note "Contributions are Acknowledged"

    Thank you for reading this page! If you need additional incentive to contribute to quacc, we plan to eventually submit a manuscript based on this code. If you contribute to quacc in a significant way (e.g. adding a new recipe or module), you will be contacted to opt-in for authorship once a draft is written. Regardless, all contributors are acknowledged in the [About](../about/contributors.md) section of the documentation.

## Steps to Contribute

To contribute to quacc, we recommend doing the following:

- [Fork the repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo) to your personal GitHub account.

- [Clone this forked repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) to your local machine, e.g. via `git clone <repo url>.git`.

- In the newly downloaded `quacc` base directory, run `pip install -e '.[dev]'` to install quacc in editable mode and with the development dependencies. For reproducibility purposes, we strongly recommend installing quacc in a fresh virtual environment.

- [Commit your changes](https://github.com/git-guides/git-commit) and [push them](https://github.com/git-guides/git-push) to your personal forked repository _in a new branch_.

- Create a [pull request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to merge your changes into the main quacc repository.

## Guidelines

### General

Please abide by the following guidelines when contributing code to quacc:

- All changes should have associated unit tests that pass. Code coverage should be maintained.

- All code should include type hints and have internally consistent documentation for the inputs and outputs.

- Use absolute paths and avoid using `os.chdir` or `os.getcwd()`. This helps ensure that quacc can be run in multi-threading mode.

- `gzip` large test files to save space.

### Writing New Recipes

!!! Tip

    In general, the easiest way to develop a new recipe is to first make one for a cheap-to-run calculator (e.g. EMT for solids, LJ for molecules), which you can then basically copy and paste for your more expensive calculator of choice.

If contributing new recipes, please abide by the following guidelines:

- Please mimic one of the other recipes as a general template for internal consistency.

- Decorate your individual compute jobs with a `#!Python @job` decorator and your flows with a `#!Python @flow` decorator. Name them accordingly.

- In general, your tests for new recipes should use a small molecule/solid (and cheap settings, where possible) to ensure that the unit tests run quickly.

### Style

In general, please try to keep the code style consistent when possible, particularly:

1. Use NumPy-style docstrings.

2. Address any relevant issues raised by the GitHub Actions test suite.

3. All Python code should be formatted with [ruff](https://github.com/astral-sh/ruff) (`ruff check --fix && ruff format`).

## Unit Tests

### General

All changes you make to quacc should be accompanied by unit tests and should not break existing tests. The full test suite will run automatically via GitHub Actions. If you are fixing a bug, include a regression test to make sure the fix you submit continues to work.

### Running Locally

To run the core test suite locally, run `pytest tests/core`.

If you wish to test the behavior of a given workflow engine, run `pytest tests/WORKFLOW_ENGINE` where `WORKFLOW_ENGINE` is the name of your workflow engine.

### Coverage

Each PR will report the code coverage once your tests pass, but if you'd like to generate a coverage report locally, you can use [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/), such as by doing `pytest --cov=quacc .`

### Monkeypatching

For recipes that involve a calculation that cannot be easily run in the test suite (e.g. the underlying calculator cannot be installed via `pip` or `conda`), you will need to make a `conftest.py` file in the recipe directory that monkeypatches the relevant functions (refer to the Gaussian recipe tests for an example). Ideally, the tests should be written in such a way that they would also pass if the executable were actually available and the `conftest.py` file were not applied. To run the test on your own machine without the `conftest.py` file applied, you can use the `--noconftest` flag when calling `pytest`. If this is not practical (e.g. because it would take too long), a separate minimal test suite can be made specifically for running on an HPC machine.

### Running a PR on an HPC Machine

It is possible to trigger the test suite on a Princeton-hosted HPC machine where the licensed executables are run without monkeypatching. If you are a member of the [@Quantum-Accelerators](https://github.com/Quantum-Accelerators), this will happen automatically. If you are not, then an admin will need to give permission. Permission is only needed if your PR is modifying a calculator or recipe.
