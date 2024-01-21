import os
from pathlib import Path

try:
    import covalent as ct
except ImportError:
    ct = None

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"


def pytest_sessionstart():
    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
    os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
    os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)

    if ct:
        ct.set_config("executors.dask.create_unique_workdir", True)
        ct.set_config("executors.local.create_unique_workdir", True)


def pytest_sessionfinish(exitstatus):
    if exitstatus == 0:
        from shutil import rmtree

        rmtree(TEST_RESULTS_DIR, ignore_errors=True)
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
