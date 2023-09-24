import os
from pathlib import Path
from shutil import rmtree

import pytest

FILE_DIR = Path(__file__).resolve().parent
TEST_RESULTS_DIR = FILE_DIR / ".test_results"
TEST_SCRATCH_DIR = FILE_DIR / ".test_scratch"


@pytest.fixture
def default_settings():
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.copy()
    DEFAULT_SETTINGS.RESULTS_DIR = TEST_RESULTS_DIR
    DEFAULT_SETTINGS.SCRATCH_DIR = TEST_SCRATCH_DIR
    return DEFAULT_SETTINGS


def pytest_sessionstart():
    os.makedirs(TEST_RESULTS_DIR, exist_ok=True)
    os.makedirs(TEST_SCRATCH_DIR, exist_ok=True)


def pytest_sessionfinish():
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
