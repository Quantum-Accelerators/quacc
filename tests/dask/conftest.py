from __future__ import annotations

import contextlib
from importlib.util import find_spec
from pathlib import Path
from shutil import rmtree

import pytest

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"

has_distributed = bool(find_spec("dask.distributed"))

if has_distributed:

    def pytest_sessionstart():

        from dask.distributed import Client, get_client
        monkeypatch = pytest.MonkeyPatch()

        file_dir = Path(__file__).parent
        monkeypatch.setenv("QUACC_CONFIG_FILE", str(file_dir / "quacc.yaml"))
        monkeypatch.setenv("QUACC_RESULTS_DIR", str(TEST_RESULTS_DIR))
        monkeypatch.setenv("QUACC_SCRATCH_DIR", str(TEST_SCRATCH_DIR))

        try:
            get_client()
        except ValueError:
            Client()

    def pytest_sessionfinish(exitstatus):
        rmtree(TEST_RESULTS_DIR, ignore_errors=True)

        if exitstatus == 0:
            from dask.distributed import default_client

            with contextlib.suppress(Exception):
                default_client().close()

            rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
