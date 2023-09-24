import importlib

import pytest
from parsl.errors import ConfigurationError


@pytest.fixture(autouse=True)
def reload_quacc(default_settings):
    importlib.reload(importlib.import_module("quacc"))
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = "parsl"
    SETTINGS.SCRATCH_DIR = default_settings.SCRATCH_DIR
    SETTINGS.RESULTS_DIR = default_settings.RESULTS_DIR


@pytest.fixture(autouse=True)
def start_server():
    from parsl.dataflow.dflow import DataFlowKernelLoader

    try:
        DataFlowKernelLoader.dfk()
    except ConfigurationError:
        import parsl

        parsl.load()
