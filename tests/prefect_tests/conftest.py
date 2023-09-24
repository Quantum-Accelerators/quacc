import importlib

import pytest

import quacc
from quacc import SETTINGS

try:
    from prefect.testing.utilities import prefect_test_harness
except ImportError:
    prefect_test_harness = None

DEFAULT_SETTINGS = SETTINGS.copy()


@pytest.mark.skipif(prefect_test_harness is None, reason="Need prefect test suite")
@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


@pytest.fixture(autouse=True)
def setup_fixture():
    importlib.reload(quacc)
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "prefect"
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR


@pytest.fixture(autouse=True)
def teardown_fixture():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR
