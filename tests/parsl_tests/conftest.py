import contextlib
import importlib

import pytest

import quacc
from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()


@pytest.fixture(autouse=True)
def setup_fixture():
    importlib.reload(quacc)
    import parsl

    from quacc import SETTINGS

    with contextlib.suppress(Exception):
        parsl.load()

    SETTINGS.WORKFLOW_ENGINE = "parsl"
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR


@pytest.fixture(autouse=True)
def teardown_fixture():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR
