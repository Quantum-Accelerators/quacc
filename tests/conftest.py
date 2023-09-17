import os
from pathlib import Path
from shutil import rmtree

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
    SETTINGS.PRIMARY_STORE = {
        "@module": "maggma.stores.mongolike",
        "@class": "MemoryStore",
        "@version": "0.56.0",
        "collection_name": "memory_db",
    }
    if not os.path.exists(SETTINGS.RESULTS_DIR):
        os.mkdir(SETTINGS.RESULTS_DIR)
    if not os.path.exists(SETTINGS.SCRATCH_DIR):
        os.mkdir(SETTINGS.SCRATCH_DIR)

    if parsl:
        try:
            parsl.load()
        except RuntimeError:
            pass


def pytest_sessionfinish():
    if os.path.exists(TEST_RESULTS_DIR):
        rmtree(TEST_RESULTS_DIR)
    if os.path.exists(TEST_SCRATCH_DIR):
        rmtree(TEST_SCRATCH_DIR)
