from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from shutil import rmtree

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"

has_ray = bool(find_spec("ray"))

if has_ray:

    def pytest_sessionstart():
        import os

        import ray

        file_dir = Path(__file__).parent
        os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
        os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
        os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)

        if not ray.is_initialized():
            ray.init(
                num_cpus=2,
                include_dashboard=False,
                ignore_reinit_error=True,
                log_to_driver=False,
                configure_logging=False,
            )

    def pytest_sessionfinish(exitstatus):
        import ray

        rmtree(TEST_RESULTS_DIR, ignore_errors=True)

        if exitstatus == 0:
            if ray.is_initialized():
                ray.shutdown()
            rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
