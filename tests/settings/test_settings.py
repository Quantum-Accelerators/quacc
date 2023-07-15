import os

from ase.build import bulk
from maggma.stores import MontyStore

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job

DEFAULT_RESULTS_STORE = SETTINGS.RESULTS_STORE
DEFAULT_GZIP_FILES = SETTINGS.GZIP_FILES
STORE = MontyStore("quacc_test_settings", database_path=".")


def teardown_function():
    SETTINGS.RESULTS_STORE = DEFAULT_RESULTS_STORE
    SETTINGS.GZIP_FILES = DEFAULT_GZIP_FILES


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
