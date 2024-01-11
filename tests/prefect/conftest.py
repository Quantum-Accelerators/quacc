import os
from pathlib import Path

import pytest

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"

try:
    import prefect
except ImportError:
    prefect = None

if prefect:

    def pytest_sessionstart():
        from prefect.testing.utilities import prefect_test_harness

        file_dir = Path(__file__).parent
        os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
        os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
        os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)

        @pytest.fixture(autouse=True, scope="module")
        def prefect_test_fixture():
            with prefect_test_harness():
                yield

    def pytest_sessionfinish():
        from shutil import rmtree

        rmtree(TEST_RESULTS_DIR, ignore_errors=True)
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
