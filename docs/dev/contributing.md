# Contributing Guide

We welcome all contributions, big or small (provided they come with unit tests!). Here are some general guidelines and tips. Note that it is never too early to open a PR if you're looking for input or feedback!

!!! Note "Contributions are Acknowledged"

    Thank you for reading this page! If you need additional incentive to contribute to quacc, we plan to eventually submit a manuscript based on this code. If you contribute to quacc in a significant way (e.g. adding a new recipe or module), you will be contacted to opt-in for authorship once a draft is written. Regardless, all contributors are acknowledged in the [About](../about/contributors.md) section of the documentation.

## Steps to Contribute

To contribute to quacc, we recommend doing the following:

- [Fork the repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo) to your personal GitHub account.

- [Clone this forked repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) to your local machine, e.g. via `git clone <repo url>.git`.

- In the newly downloaded `quacc` base directory, run `pip install -e .[dev]` to install quacc in editable mode and with the development dependencies. For reproducibility purposes, we strongly recommend installing quacc in a fresh virtual environment.

- [Commit your changes](https://github.com/git-guides/git-commit) and [push them](https://github.com/git-guides/git-push) to your personal forked repository _in a new branch_.

- Create a [pull request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to merge your changes into the main quacc repository.

## Guidelines

### General

Please abide by the following guidelines when contributing code to quacc:

- All changes should have associated unit tests that pass. Code coverage should be maintained.

- All code should include type hints and have internally consistent documentation for the inputs and outputs.

- Use absolute paths where possible and avoid using `os.chdir`. This helps ensure that quacc can be run in multi-threading mode.

- `gzip` large test files to save space.

### Writing New Recipes

!!! Tip

    In general, the easiest way to develop a new recipe is to first make one for a cheap-to-run calculator (e.g. EMT for solids, LJ for molecules), which you can then basically copy and paste for your more expensive calculator of choice.

If contributing new recipes, please abide by the following guidelines:

- Please mimic one of the other recipes as a general template for internal consistency.

- Decorate your individual compute jobs with a `#!Python @job` decorator and your flows with a `#!Python @flow` decorator. Name them accordingly.

- Add your new recipe to the [List of Recipes](../user/recipes/recipes_list.md).

- In general, your tests for new recipes should use a small molecule/solid (and cheap settings, where possible) to ensure that the unit tests run quickly.

### Style

In general, please try to keep the code style consistent when possible, particularly:

1. Use NumPy-style docstrings.

2. Address any relevant issues raised by the GitHub Actions test suite.

3. All Python code should be formatted with [isort](https://github.com/PyCQA/isort) (`isort .`), [black](https://github.com/psf/black) (`black .`), and [ruff](https://github.com/astral-sh/ruff) (`ruff . --fix`), although this will be corrected automatically when merged.

## Unit Tests

### General

All changes you make to quacc should be accompanied by unit tests and should not break existing tests. The full test suite will run automatically via GitHub Actions. If you are fixing a bug, include a regression test to make sure the fix you submit continues to work.

### Running Locally

To run the core test suite locally, run `pytest tests/local`.

If you wish to test the behavior of a given workflow engine, run `pytest tests/WORKFLOW_ENGINE` where `WORKFLOW_ENGINE` is the name of your workflow engine.

### Coverage

Each PR will report the code coverage once your tests pass, but if you'd like to generate a coverage report locally, you can use [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/), such as by doing `pytest --cov=quacc .`.

!!! Note

    The codecov check on PRs updates as each unit test finishes, so don't be surprised to see an "X" if tests are still running.

### Monkeypatching

For recipes that involve a calculation that cannot be easily run in the test suite (e.g. the underlying calculator cannot be installed via `pip` or `conda`), you will need to make a `conftest.py` file in the recipe directory that monkeypatches the relevant functions (refer to the Gaussian recipe tests for an example). The tests should be written in such a way that they would also pass if the executable were actually available and the `conftest.py` file were not applied. To run the test on your own machine without the `conftest.py` file applied, you can use the `--noconftest` flag when calling `pytest`.

### Jenkins

When you submit a PR, a job is submitted on the [Adroit cluster](https://researchcomputing.princeton.edu/systems/adroit) at Princeton that will automatically run the test suite on HPC resources without relying on monkeypatching via `conftest.py` files (i.e. using the licensed executables). This is done via a [Jenkins server](https://www.jenkins.io/). If you are adding recipes for a calculator that is not yet part of quacc and cannot be installed via `pip` or `conda`, we will add the necessary executables to the [Jenkins build pipeline](https://github.com/Quantum-Accelerators/quacc/tree/main/tests/jenkins) so that production-quality tests will be run automatically (provided we have a license for it). This is an important reason for why the unit tests should be done with a small molecule/solid.
