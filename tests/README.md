This is the test suite for quacc, which is executed via the `pytest` utility.

Note that you aren't meant to run `pytest` from the root of the `tests` directory. Instead, you need to run `pytest` from within one of the subdirectories found in `tests`. This is because each subdirectory contains a `conftest.py` file that is used to set up a given test environment at runtime, and they will otherwise clash with each other.

For instance, to run the core test suite:

```bash
pytest tests/local # or `pytest .` from within `tests/local`
```

Similarly, to run the test suite for Parsl-specific functionality:

```bash
pytest tests/parsl # or `pytest .` from within `tests/parsl`
```
