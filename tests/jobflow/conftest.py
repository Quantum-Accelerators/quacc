from __future__ import annotations

from pathlib import Path
from shutil import rmtree

import pytest

TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"


@pytest.fixture
def jobflow_output():
    """Return the latest stored output for a Jobflow job."""

    def get_output(responses, job):
        indexed_responses = responses[job.uuid]
        return indexed_responses[max(indexed_responses)].output

    return get_output


def pytest_sessionstart():
    import os

    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
    os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)
    os.environ["JOBFLOW_CONFIG_FILE"] = str(file_dir / "jobflow.yaml")


def pytest_sessionfinish(exitstatus):
    if exitstatus == 0:
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
