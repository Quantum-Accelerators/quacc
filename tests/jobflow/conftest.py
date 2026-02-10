from __future__ import annotations

from pathlib import Path
from shutil import rmtree

TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"


def pytest_sessionstart():
    import os

    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
    os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)
    os.environ["JOBFLOW_CONFIG_FILE"] = str(file_dir / "jobflow.yaml")


def pytest_sessionfinish(exitstatus):
    if exitstatus == 0:
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)


import pytest


@pytest.fixture(autouse=True)
def reset_contextvars():
    from quacc.wflow_tools.context import _directory_context, _execution_context

    token_exec = _execution_context.set(())
    token_dir = _directory_context.set("")

    yield

    # Restore whatever was there before the test
    _execution_context.reset(token_exec)
    _directory_context.reset(token_dir)
