import os
from shutil import rmtree

import pytest
from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job

try:
    import montydb
except:
    montydb = None
DEFAULT_SETTINGS = SETTINGS.copy()


def teardown_function():
    SETTINGS.RESULTS_STORE = DEFAULT_SETTINGS.RESULTS_STORE
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


@pytest.mark.skipif(montydb is None, reason="MontyDB not installed")
def test_store():
    store = MontyStore("quacc_test_settings", database_path=".")
    SETTINGS.RESULTS_STORE = store.to_json()
    atoms = bulk("Cu")
    static_job(atoms)
    store.connect()
    assert store.count() == 1
    store.close()


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
