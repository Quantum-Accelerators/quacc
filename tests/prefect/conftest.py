def pytest_sessionstart():
    import os
    from pathlib import Path

    from quacc import SETTINGS

    file_dir = Path(__file__).resolve().parent
    test_results_dir = file_dir / ".test_results"
    test_scratch_dir = file_dir / ".test_scratch"

    SETTINGS.WORKFLOW_ENGINE = "prefect"
    SETTINGS.RESULTS_DIR = test_results_dir
    SETTINGS.SCRATCH_DIR = test_scratch_dir
    os.makedirs(test_results_dir, exist_ok=True)
    os.makedirs(test_scratch_dir, exist_ok=True)


def pytest_sessionfinish():
    from pathlib import Path
    from shutil import rmtree

    file_dir = Path(__file__).resolve().parent
    test_results_dir = file_dir / ".test_results"
    test_scratch_dir = file_dir / ".test_scratch"

    rmtree(test_results_dir, ignore_errors=True)
    rmtree(test_scratch_dir, ignore_errors=True)
