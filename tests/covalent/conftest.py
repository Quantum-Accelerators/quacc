from pathlib import Path

from quacc import SETTINGS


def pytest_sessionstart():
    test_results_dir = Path.cwd() / ".test_results"
    test_scratch_dir = Path.cwd() / ".test_scratch"
    SETTINGS.WORKFLOW_ENGINE = "covalent"
    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
