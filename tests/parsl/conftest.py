import os
from pathlib import Path

TEST_RESULTS_DIR = Path(__file__).parent / "test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "test_scratch"
TEST_RUNINFO = Path(__file__).parent / "runinfo"


def pytest_sessionstart():
    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / ".quacc.yaml")


def pytest_sessionfinish():
    from shutil import rmtree

    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
    rmtree(TEST_RUNINFO, ignore_errors=True)
