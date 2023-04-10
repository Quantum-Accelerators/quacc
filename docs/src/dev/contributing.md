
# Contributing

To contribute to QuAcc, we recommend installing QuAcc in editable mode so that your changes are reflected every time you run QuAcc. This can be done by using the `-e` flag when running `pip install`. Also, you will need to install the development dependencies via `pip install .[dev]`.

To ensure a consistent formatting style, we use `isort` and `black` to format the code. In the base directory of `quacc`, simply run `isort .` and then `black .` (in that order) to format the code. This will ensure your contributions pass the linting check.

To ensure that any changes you have made do not break unit tests, run `pytest` from the the `quacc/tests` directory. Any new features that you propose to add to QuAcc must have accompanying unit tests.

If the recipes you're adding are based on a code that can be readily installed via `pip` or `conda` (e.g. tblite, xtb-python, DFTB+, Psi4), then you can run the unit tests with these codes directly. Preferably, you should use a small molecule or solid so the unit tests run quickly. If the recipes you're adding are proprietary or not available via `pip` or `conda` (e.g. Gaussian, GULP), then you will need to monkeypatch certain functions to change their behavior during testing. For instance, we do not want to run VASP directly during unit tests and have mocked the `atoms.get_potential_energy()` function to always return a value of -1.0 during unit tests. Any mocked functions can be found in the `conftest.py` files of the testing directory.

In short, when contributing new recipes and unit tests, please follow the existing style and structure of the code where possible.
