from pathlib import Path

from maggma.stores import MemoryStore

from quacc import SETTINGS


def pytest_sessionstart():
    test_results_dir = Path.cwd() / ".test_results"
    test_scratch_dir = Path.cwd() / ".test_scratch"
    SETTINGS.PRIMARY_STORE = MemoryStore()
    SETTINGS.WORKFLOW_ENGINE = "local"
    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
