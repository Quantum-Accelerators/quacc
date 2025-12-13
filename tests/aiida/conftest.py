from __future__ import annotations

import subprocess
from importlib.util import find_spec
from pathlib import Path
from shutil import rmtree

TEST_RESULTS_DIR = Path(__file__).parent / "_test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "_test_scratch"

has_aiida = bool(find_spec("aiida_workgraph"))

if has_aiida:
    subprocess.run(
        ["verdi", "presto", "--profile-name", "test_profile"],
        capture_output=True,
        text=True,
        check=False,
    )

    def pytest_sessionstart():
        import os

        from aiida import load_profile

        file_dir = Path(__file__).parent
        os.environ["QUACC_CONFIG_FILE"] = str(file_dir / "quacc.yaml")
        os.environ["QUACC_RESULTS_DIR"] = str(TEST_RESULTS_DIR)
        os.environ["QUACC_SCRATCH_DIR"] = str(TEST_SCRATCH_DIR)

        load_profile("test_profile", allow_switch=True)

    def pytest_sessionfinish(exitstatus):
        rmtree(TEST_RESULTS_DIR, ignore_errors=True)
        if exitstatus == 0:
            rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
