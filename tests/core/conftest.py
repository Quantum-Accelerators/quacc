from __future__ import annotations

import os
from pathlib import Path
from shutil import rmtree

TEST_RENDERS_DIR = Path(__file__).parent / "_test_renders"
TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"


def pytest_sessionstart():
    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / ".quacc.yaml")
    os.environ["QUACC_RENDERS_DIR"] = str(TEST_RENDERS_DIR)
    os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
    os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)


def pytest_sessionfinish(exitstatus):
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    if exitstatus == 0:
        rmtree(TEST_RENDERS_DIR, ignore_errors=True)
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
