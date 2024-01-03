import os
from pathlib import Path
from shutil import rmtree

TEST_RESULTS_DIR = Path(__file__).parent / "test_results"
TEST_SCRATCH_DIR = Path(__file__).parent / "test_scratch"

try:
    import dask.distributed

    has_import = True
except ImportError:
    has_import = False

if has_import:

    def pytest_sessionstart():
        from dask.distributed import Client, default_client

        file_dir = Path(__file__).parent
        os.environ["QUACC_CONFIG_FILE"] = str(file_dir / ".quacc.yaml")
        try:
            default_client()
        except ValueError:
            Client()

    def pytest_sessionfinish():
        from dask.distributed import default_client

        try:
            default_client().close()
        except Exception:
            pass

        rmtree(TEST_RESULTS_DIR, ignore_errors=True)
        rmtree(TEST_SCRATCH_DIR, ignore_errors=True)
