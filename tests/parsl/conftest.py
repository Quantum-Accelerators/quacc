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
    from parsl.dataflow.memoization import BasicMemoizer


def pytest_sessionstart():
    import os

    if has_parsl:
        parsl.load(
            Config(
                dependency_resolver=DEEP_DEPENDENCY_RESOLVER,
                run_dir=str(TEST_RUNINFO),
                memoizer=BasicMemoizer(checkpoint_mode="task_exit"),
            )
        )
    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
    os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
    os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)


def pytest_sessionfinish(exitstatus):
    if has_parsl:
        parsl.clear()
    rmtree(TEST_RESULTS_DIR, ignore_errors=True)
    if exitstatus == 0:
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
        rmtree(TEST_RUNINFO, ignore_errors=True)
