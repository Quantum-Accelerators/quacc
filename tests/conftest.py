import os
from pathlib import Path
from shutil import rmtree

from quacc import SETTINGS

FILE_DIR = Path(__file__).resolve().parent

test_results_dir = os.path.join(FILE_DIR, ".results_dir")
test_scratch_dir = os.path.join(FILE_DIR, ".scratch_dir")


def pytest_sessionstart():
    os.makedirs(test_results_dir, exist_ok=True)
    os.makedirs(test_scratch_dir, exist_ok=True)
    SETTINGS.WORKFLOW_ENGINE = None
    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir


def pytest_sessionfinish():
    if os.path.exists(test_results_dir):
        rmtree(test_results_dir)
    if os.path.exists(test_scratch_dir):
        rmtree(test_scratch_dir)
