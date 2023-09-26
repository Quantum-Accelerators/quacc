import os
from pathlib import Path
from shutil import rmtree

FILE_DIR = Path(__file__).resolve().parent
TEST_RESULTS_DIR = FILE_DIR / ".test_results"
TEST_SCRATCH_DIR = FILE_DIR / ".test_scratch"


def pytest_sessionstart():
    from quacc import SETTINGS

    SETTINGS.RESULTS_DIR = TEST_RESULTS_DIR
    SETTINGS.SCRATCH_DIR = TEST_SCRATCH_DIR


def pytest_sessionfinish():
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
