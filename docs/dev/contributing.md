# Contributing Guide

We welcome all contributions, big or small (provided they come with unit tests!). Here are some general guidelines and tips. Note that it is never too early to open a PR if you're looking for input or feedback!

!!! Note "Contributions are Acknowledged"

    Thank you for reading this page! If you need additional incentive to contribute to quacc, we plan to eventually submit a manuscript based on this code. If you contribute to quacc in a significant way (e.g. adding a new recipe or module), you will be contacted to opt-in for authorship once a draft is written. Regardless, all contributors are acknowledged in the [About](../about/contributors.md) section of the documentation.

## General Tips

When developing new recipes, it is often helpful to start from an existing example. In general, we recommend referring to `quacc.recipes.emt` as a starting point.

## Steps to Contribute

To contribute to quacc, we recommend doing the following:

- [Fork the repository](https://docs.github.com/en/get-started/quickstart/fork-a-repo) to your personal GitHub account.

- [Clone this forked repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) to your local machine, e.g. via `git clone <repo url>.git`.

- In the newly downloaded `quacc` base directory, run `pip install -e .[dev]` to install quacc in editable mode and with the development dependencies. If you need to strictly match the versions of dependencies used in the GitHub actions test suite, you can find `requirements.txt` files in the `tests` directory.

- [Commit your changes](https://github.com/git-guides/git-commit) and [push them](https://github.com/git-guides/git-push) to your personal forked repository _in a new branch_.

- Create a [pull request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to merge your changes into the main quacc repository.

For reproducibility purposes, we strongly recommend installing quacc in a fresh virtual environment.

## Guidelines

### General

Please abide by the following guidelines when contributing code to quacc:

- All changes should have associated unit tests that pass. Code coverage should be maintained.

- All code should include type hints and have internally consistent documentation for the inputs and outputs.

- Only define multi-step workflows if they go beyond simply stitching together existing functions or if they are widely used in other recipes. Otherwise, just define the individual functions.

- Use absolute paths where possible and avoid using `os.chdir`. This helps ensure that quacc can be run in a multi-threading mode.

- Put your imports inside each function in your test suite

- `gzip` large test files to save space.

### New Recipes

If contributing new recipes, please abide by the following guidelines:

- Please mimic one of the other recipes as a general template for internal consistency.

- Decorate your individual compute jobs with a `#!Python @job` decorator or your flow with a `#!Python @flow` decorator and name them accordingly.

In general, the easiest way to develop a new recipe is to first make one for a cheap-to-run calculator (e.g. EMT for solids, LJ for molecules), which you can then basically copy and paste for your more expensive calculator of choice.

### Style

In general, please try to keep the code style consistent when possible, particularly:

1. Use NumPy-style Docstrings.

2. Address any relevant style issues raised by the "Deepsource: Python," if present.

3. All Python code should be formatted with [isort](https://github.com/PyCQA/isort), [black](https://github.com/psf/black), and [ruff](https://github.com/astral-sh/ruff), although this will be corrected automatically when merged.

## Unit Tests

### General

All changes you make to quacc should be accompanied by unit tests and should not break existing tests. The full test suite will run on automatically via GitHub actions.

### Running Locally

To run the core test suite locally, run `pytest tests/local` (which is the same as doing `pytest` without any extra arguments).

If you wish to test the behavior of a given workflow engine, run `pytest tests/WORKFLOW_ENGINE` where `WORKFLOW_ENGINE` is the name of your workflow engine.

### Coverage

Each PR will report the coverage once your tests pass, but if you'd like to generate a coverage report locally, you can use [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/), such as by doing `pytest --cov=quacc .`.

!!! Note

    The codecov check on pull requests update as each unit test finishes, so don't be surprised to see an "X" if tests are still running.

### Adding New Recipes

If you are adding recipes based on a code that can be readily installed via `pip` or `conda` (e.g. tblite, DFTB+, Psi4), then you can run these codes directly in the test suite. Preferably, you should use a small molecule or solid and cheap method so the unit tests run quickly.

If the recipes you're adding are proprietary or not available via `pip` or `conda` (e.g. Gaussian, GULP), then you will need to [monkeypatch](https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html) certain functions to change their behavior during testing.
