import os

import pytest
from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job

try:
    import montydb
except:
    montydb = None
DEFAULT_SETTINGS = SETTINGS.copy()
STORE = MontyStore("quacc_test_settings", database_path=".")


def teardown_function():
    SETTINGS.RESULTS_STORE = DEFAULT_SETTINGS.RESULTS_STORE
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR


@pytest.skipif(montydb is None, reason="MontyDB not installed")
def test_store():
    SETTINGS.RESULTS_STORE = STORE.to_json()
    atoms = bulk("Cu")
    static_job(atoms)
    STORE.connect()
    assert STORE.count() == 1
    STORE.close()


def test_results_dir():
    atoms = bulk("Cu")
    relax_job(atoms)
    assert "opt.traj.gz" in os.listdir(os.getcwd())
    os.remove("opt.traj.gz")
    SETTINGS.GZIP_FILES = False
    relax_job(atoms)
    assert "opt.traj" in os.listdir(os.getcwd())
    os.remove("opt.traj")


def test_create_unique_workdir():
    atoms = bulk("Cu")
    relax_job(atoms)
    assert "quacc_" not in os.listdir(os.getcwd())
    SETTINGS.CREATE_UNIQUE_WORKDIR = True
    relax_job(atoms)
    assert "quacc_" in os.listdir(os.getcwd())
