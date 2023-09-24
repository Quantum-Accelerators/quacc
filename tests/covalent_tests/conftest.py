import importlib
import subprocess

import pytest


@pytest.fixture(autouse=True)
def reload_quacc(default_settings):
    importlib.reload(importlib.import_module("quacc"))
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "covalent"
    SETTINGS.SCRATCH_DIR = default_settings.SCRATCH_DIR
    SETTINGS.RESULTS_DIR = default_settings.RESULTS_DIR


@pytest.fixture(autouse=True)
def start_server():
    subprocess.run(["covalent", "start"], check=True)
