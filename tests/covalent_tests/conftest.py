import importlib
import subprocess

import pytest
import requests


@pytest.fixture(autouse=True)
def reload_quacc(default_settings):
    importlib.reload(importlib.import_module("quacc"))
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "covalent"
    SETTINGS.SCRATCH_DIR = default_settings.SCRATCH_DIR
    SETTINGS.RESULTS_DIR = default_settings.RESULTS_DIR


@pytest.fixture(autouse=True)
def start_server():
    import covalent as ct

    address = f"http://{ct.get_config('dispatcher.address')}:{str(ct.get_config('dispatcher.port'))}"
    try:
        requests.session().get(address)
    except requests.exceptions.ConnectionError:
        subprocess.run(["covalent", "start"], check=True)
