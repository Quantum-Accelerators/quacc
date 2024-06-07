from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from shutil import rmtree

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"
TEST_RUNINFO = Path(__file__).parent / "runinfo"

has_parsl = bool(find_spec("parsl"))
if has_parsl:
    import parsl
    from parsl.config import Config
    from parsl.dataflow.dependency_resolvers import DEEP_DEPENDENCY_RESOLVER


def pytest_sessionstart():
    import pytest
    monkeypatch = pytest.MonkeyPatch()

    if parsl:
        parsl.load(Config(dependency_resolver=DEEP_DEPENDENCY_RESOLVER))
    file_dir = Path(__file__).parent
    monkeypatch.setenv("QUACC_CONFIG_FILE", str(file_dir / "quacc.yaml"))
    monkeypatch.setenv("QUACC_RESULTS_DIR", str(TEST_RESULTS_DIR))
    monkeypatch.setenv("QUACC_SCRATCH_DIR", str(TEST_SCRATCH_DIR))


def pytest_sessionfinish(exitstatus):
    if parsl:
        parsl.clear()
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    if exitstatus == 0:
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
        rmtree(TEST_RUNINFO, ignore_errors=True)
