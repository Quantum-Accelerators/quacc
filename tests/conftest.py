import os
from pathlib import Path
from shutil import rmtree

from quacc import SETTINGS

FILE_DIR = Path(__file__).resolve().parent

DEFAULT_SETTINGS = SETTINGS.copy()

test_results_dir = os.path.join(FILE_DIR, ".results_dir")
test_scratch_dir = os.path.join(FILE_DIR, ".scratch_dir")

WFLOW_ENGINE = SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None


def pytest_sessionstart():
    os.makedirs(test_results_dir, exist_ok=True)
    os.makedirs(test_scratch_dir, exist_ok=True)
    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
    if WFLOW_ENGINE == "parsl":
        import parsl

        parsl.load()


def pytest_sessionfinish():
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR
    rmtree(test_results_dir)
    rmtree(test_scratch_dir)
