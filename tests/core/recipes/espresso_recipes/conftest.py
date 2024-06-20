from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def ESPRESSO_PARALLEL_CMD():
    from shutil import which

    import psutil

    from quacc import change_settings

    if which("mpirun") and psutil.cpu_count(logical=False) >= 2:
        change_settings({"ESPRESSO_PARALLEL_CMD": "mpirun -np 2"})
    else:
        change_settings({"ESPRESSO_PARALLEL_CMD": ""})
