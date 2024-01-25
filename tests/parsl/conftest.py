import os
from pathlib import Path

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"
TEST_RUNINFO = Path(__file__).parent / "runinfo"

try:
    import parsl
except ImportError:
    parsl = None


def pytest_sessionstart():
    if parsl:
        from parsl.channels import LocalChannel
        from parsl.config import Config
        from parsl.executors import HighThroughputExecutor
        from parsl.launchers import SimpleLauncher
        from parsl.providers import LocalProvider

        config = Config(
            executors=[
                HighThroughputExecutor(
                    label="test_htex_local",
                    worker_debug=True,
                    cores_per_worker=1,
                    provider=LocalProvider(
                        channel=LocalChannel(),
                        init_blocks=1,
                        max_blocks=1,
                        launcher=SimpleLauncher(),
                    ),
                )
            ],
            strategy="none",
        )
        parsl.load(config)
    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
    os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
    os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)


def pytest_sessionfinish(exitstatus):
    if parsl:
        parsl.clear()
    if exitstatus == 0:
        from shutil import rmtree

        rmtree(TEST_RESULTS_DIR, ignore_errors=True)
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
        rmtree(TEST_RUNINFO, ignore_errors=True)
