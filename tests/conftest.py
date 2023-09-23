import contextlib
import os
from pathlib import Path
from shutil import rmtree

from maggma.stores import MemoryStore

try:
    import parsl
except ImportError:
    parsl = None
FILE_DIR = Path(__file__).resolve().parent
TEST_RESULTS_DIR = FILE_DIR / ".test_results"
TEST_SCRATCH_DIR = FILE_DIR / ".test_scratch"


def pytest_sessionstart():
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "local"
    SETTINGS.RESULTS_DIR = str(TEST_RESULTS_DIR)
    SETTINGS.SCRATCH_DIR = str(TEST_SCRATCH_DIR)
    SETTINGS.PRIMARY_STORE = MemoryStore().to_json()
    os.makedirs(SETTINGS.RESULTS_DIR, exist_ok=True)
    os.makedirs(SETTINGS.SCRATCH_DIR, exist_ok=True)

    if parsl:
        with contextlib.suppress(RuntimeError):
            parsl.load()


def pytest_sessionfinish():
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
