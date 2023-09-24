import contextlib
import os
import subprocess
from pathlib import Path
from shutil import rmtree

from quacc import SETTINGS

try:
    import covalent as ct
except ImportError:
    ct = None
try:
    import parsl
except ImportError:
    parsl = None
FILE_DIR = Path(__file__).resolve().parent
TEST_RESULTS_DIR = FILE_DIR / ".test_results"
TEST_SCRATCH_DIR = FILE_DIR / ".test_scratch"


def pytest_sessionstart():
    SETTINGS.RESULTS_DIR = str(TEST_RESULTS_DIR)
    SETTINGS.SCRATCH_DIR = str(TEST_SCRATCH_DIR)

    os.makedirs(SETTINGS.RESULTS_DIR, exist_ok=True)
    os.makedirs(SETTINGS.SCRATCH_DIR, exist_ok=True)

    if ct:
        subprocess.run(["covalent", "start"], check=True)
    if parsl:
        with contextlib.suppress(Exception):
            parsl.load()


def pytest_sessionfinish():
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
