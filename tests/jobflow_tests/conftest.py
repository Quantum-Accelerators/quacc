import importlib

import pytest

import quacc
from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()


@pytest.fixture(autouse=True)
def setup_fixture():
    importlib.reload(quacc)
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "jobflow"
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR


@pytest.fixture(autouse=True)
def teardown_fixture():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR
