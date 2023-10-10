import pytest


@pytest.fixture
def default_settings():
    from quacc import SETTINGS

    return SETTINGS.copy()


def pytest_sessionstart(test_results_dir, test_scratch_dir):
    import os
    from pathlib import Path

    from quacc import SETTINGS

    file_dir = Path(__file__).resolve().parent
    test_results_dir = file_dir / ".test_results"
    test_scratch_dir = file_dir / ".test_scratch"

    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
    os.makedirs(test_results_dir, exist_ok=True)
    os.makedirs(test_scratch_dir, exist_ok=True)


def pytest_sessionfinish(test_results_dir, test_scratch_dir):
    from pathlib import Path
    from shutil import rmtree

    file_dir = Path(__file__).resolve().parent
    test_results_dir = file_dir / ".test_results"
    test_scratch_dir = file_dir / ".test_scratch"

    rmtree(test_results_dir, ignore_errors=True)
    rmtree(test_scratch_dir, ignore_errors=True)
