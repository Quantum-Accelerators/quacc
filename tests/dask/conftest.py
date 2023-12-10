import os
from pathlib import Path


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
    except:
        pass
