import importlib

import pytest


@pytest.fixture(autouse=True, scope="module")
def reload_jobflow_quacc(default_settings):
    importlib.reload(importlib.import_module("quacc"))
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "jobflow"
    SETTINGS.SCRATCH_DIR = default_settings.SCRATCH_DIR
    SETTINGS.RESULTS_DIR = default_settings.RESULTS_DIR
