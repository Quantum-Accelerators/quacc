import os
from shutil import rmtree

from ase.build import bulk
from maggma.stores import MemoryStore

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job

DEFAULT_SETTINGS = SETTINGS.copy()


def teardown_function():
    SETTINGS.RESULTS_STORE = DEFAULT_SETTINGS.RESULTS_STORE
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or "quacc_" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_store():
    store = MemoryStore()
    SETTINGS.RESULTS_STORE = store.to_json()
    atoms = bulk("Cu")
    static_job(atoms)
    with store:
        assert store.count() == 1


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
    assert "quacc_" not in os.listdir(os.getcwd())[1]
    SETTINGS.CREATE_UNIQUE_WORKDIR = True
    relax_job(atoms)
    assert "quacc_" in os.listdir(os.getcwd())[1]
