import pytest


@pytest.fixture()
def test_results_dir():
    from pathlib import Path

    FILE_DIR = Path(__file__).resolve().parent
    return FILE_DIR / ".test_results"


@pytest.fixture()
def test_scratch_dir():
    from pathlib import Path

    FILE_DIR = Path(__file__).resolve().parent

    return FILE_DIR / ".test_scratch"


@pytest.fixture
def default_settings():
    from quacc import SETTINGS

    return SETTINGS.copy()


def pytest_sessionstart():
    import os

    from quacc import SETTINGS

    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
    os.makedirs(SETTINGS.RESULTS_DIR, exist_ok=True)
    os.makedirs(SETTINGS.SCRATCH_DIR, exist_ok=True)


def pytest_sessionfinish():
    from shutil import rmtree

    rmtree(test_results_dir, ignore_errors=True)
    rmtree(test_scratch_dir, ignore_errors=True)
