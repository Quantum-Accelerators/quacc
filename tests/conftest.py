import os
from pathlib import Path
from shutil import rmtree

FILE_DIR = Path(__file__).resolve().parent
TEST_RESULTS_DIR = FILE_DIR / ".test_results"
TEST_SCRATCH_DIR = FILE_DIR / ".test_scratch"


def pytest_sessionstart():
    from quacc import SETTINGS

    SETTINGS.RESULTS_DIR = str(TEST_RESULTS_DIR)
    SETTINGS.SCRATCH_DIR = str(TEST_SCRATCH_DIR)
    if not os.path.exists(SETTINGS.RESULTS_DIR):
        os.mkdir(SETTINGS.RESULTS_DIR)
    if not os.path.exists(SETTINGS.SCRATCH_DIR):
        os.mkdir(SETTINGS.SCRATCH_DIR)

    WFLOW_ENGINE = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )

    if WFLOW_ENGINE == "parsl":
        import parsl

        parsl.load()


def pytest_sessionfinish():
    if os.path.exists(TEST_RESULTS_DIR):
        rmtree(TEST_RESULTS_DIR)
    if os.path.exists(TEST_SCRATCH_DIR):
        rmtree(TEST_SCRATCH_DIR)
