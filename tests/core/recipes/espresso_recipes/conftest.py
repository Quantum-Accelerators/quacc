import pytest


@pytest.fixture
def DEFAULT_PARALLEL_INFO():
    from shutil import which

    import psutil

    if which("mpirun") and psutil.cpu_count(logical=False) >= 2:
        return {"binary": "mpirun", "-np": 2}
    else:
        return {}
