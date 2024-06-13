# Overview

This is the test suite for quacc, which is executed via the `pytest` utility.

## Running the Tests

Note that you aren't meant to run `pytest` from the root of the `tests` directory. Instead, you need to run `pytest` from within one of the subdirectories found in `tests`. This is because each subdirectory contains a `conftest.py` file that is used to set up a given test environment at runtime, and they will otherwise clash with each other.

For instance, to run the core test suite:

```bash
pytest tests/core # or `pytest .` from within `tests/core`
```

Similarly, to run the test suite for Parsl-specific functionality:

```bash
pytest tests/parsl # or `pytest .` from within `tests/parsl`
```

## Troubleshooting Tests

If you ever need to inspect files from a failed test run, you can find them in the `_quacc_scratch` and `_quacc_results` directories that are created in the current working directory.

## Requirements

The `requirements.txt` files in this directory are those used on GitHub Actions and corresponding to the optional dependency groups listed in `pyproject.toml`.
