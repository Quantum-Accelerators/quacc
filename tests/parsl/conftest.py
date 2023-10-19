import os
from pathlib import Path

import pytest

from quacc import SETTINGS


@pytest.fixture(scope="session", autouse=True)
def set_settings():
    file_dir = Path(__file__).resolve().parent
    test_results_dir = file_dir / ".test_results"
    test_scratch_dir = file_dir / ".test_scratch"
    SETTINGS.WORKFLOW_ENGINE = "parsl"
    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
    os.makedirs(test_results_dir, exist_ok=True)
    os.makedirs(test_scratch_dir, exist_ok=True)
