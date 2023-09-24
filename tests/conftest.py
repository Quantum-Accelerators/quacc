import contextlib
import os
import subprocess
from pathlib import Path
from shutil import rmtree

from maggma.stores import MemoryStore

FILE_DIR = Path(__file__).resolve().parent
TEST_RESULTS_DIR = FILE_DIR / ".test_results"
TEST_SCRATCH_DIR = FILE_DIR / ".test_scratch"


def pytest_sessionstart():
    from quacc import SETTINGS

    SETTINGS.RESULTS_DIR = str(TEST_RESULTS_DIR)
    SETTINGS.SCRATCH_DIR = str(TEST_SCRATCH_DIR)
    SETTINGS.PRIMARY_STORE = MemoryStore().to_json()
    os.makedirs(SETTINGS.RESULTS_DIR, exist_ok=True)
    os.makedirs(SETTINGS.SCRATCH_DIR, exist_ok=True)

    if SETTINGS.WORKFLOW_ENGINE == "covalent":
        subprocess.run(["covalent", "start"])
    elif SETTINGS.WORKFLOW_ENGINE == "parsl":
        import parsl

        with contextlib.suppress(Exception):
            parsl.load()


def pytest_sessionfinish():
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
