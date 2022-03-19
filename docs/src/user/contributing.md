
# Contributing

To contribute to QuAcc, we recommend installing QuAcc in editable mode so that your changes are reflected every time you run QuAcc. This can be done by using the `-e` flag when running `pip install`. Also, you will need to install the development dependencies via `pip install .[dev]`.

To ensure a consistent formatting style, we use `isort` and `black` to format the code. In the base directory of `quacc`, simply run `isort .` and then `black .` (in that order) to format the code. This will ensure your contributions pass the linting check.

To ensure that any changes you have made do not break unit tests, run `pytest` from the the `quacc/tests` directory.

Any new features that you propose to add to QuAcc must have accompanying unit tests. In certain cases, we have monkeypatched certain functions to change their behavior during testing. For instance, we do not want to run VASP directly during unit tests and have mocked the `atoms.get_potential_energy()` function to always return a value of -1.0 during unit tests. Any mocked functions can be found in the `conftest.py` files of the testing directory.
