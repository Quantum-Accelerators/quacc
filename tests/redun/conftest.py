from __future__ import annotations

from pathlib import Path
from shutil import rmtree

import pytest

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"


def pytest_sessionstart():
    file_dir = Path(__file__).parent
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setenv("QUACC_CONFIG_FILE", str(file_dir / "quacc.yaml"))
    monkeypatch.setenv("QUACC_RESULTS_DIR", str(TEST_RESULTS_DIR))
    monkeypatch.setenv("QUACC_SCRATCH_DIR", str(TEST_SCRATCH_DIR))


def pytest_sessionfinish(exitstatus):
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    if exitstatus == 0:
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
