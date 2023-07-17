import os
from glob import glob

from ase.build import bulk
from maggma.stores import MemoryStore

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job
from quacc.settings import QuaccSettings

DEFAULT_SETTINGS = SETTINGS.copy()


def setup_function():
    SETTINGS.PRIMARY_STORE = None
    SETTINGS.GZIP_FILES = True
    SETTINGS.CREATE_UNIQUE_WORKDIR = False


def teardown_function():
    SETTINGS.PRIMARY_STORE = DEFAULT_SETTINGS.PRIMARY_STORE
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR


def test_file(monkeypatch, tmpdir):
    tmpdir.chdir()

    assert QuaccSettings().GZIP_FILES is True

    with open("quacc_test.yaml", "w") as f:
        f.write("GZIP_FILES: false")
    monkeypatch.setenv(
        "QUACC_CONFIG_FILE", os.path.join(os.getcwd(), "quacc_test.yaml")
    )

    assert QuaccSettings().GZIP_FILES is False
    os.remove("quacc_test.yaml")


def test_store():
    store = MemoryStore()
    SETTINGS.PRIMARY_STORE = store.to_json()
    atoms = bulk("Cu")
    static_job(atoms)


def test_results_dir(tmpdir):
    tmpdir.chdir()

    atoms = bulk("Cu")
    relax_job(atoms)
    assert "opt.traj.gz" in os.listdir(os.getcwd())
    os.remove("opt.traj.gz")
    SETTINGS.GZIP_FILES = False
    relax_job(atoms)
    assert "opt.traj" in os.listdir(os.getcwd())
    os.remove("opt.traj")


def test_env_var(monkeypatch):
    monkeypatch.setenv("QUACC_SCRATCH_DIR", "/my/scratch/dir")
    assert QuaccSettings().SCRATCH_DIR == "/my/scratch/dir"


def test_yaml(tmpdir, monkeypatch):
    tmpdir.chdir()

    with open("quacc_test.yaml", "w") as f:
        f.write('SCRATCH_DIR: "/my/new/scratch/dir"')
    monkeypatch.setenv("QUACC_CONFIG_FILE", "quacc_test.yaml")
    assert QuaccSettings().SCRATCH_DIR == "/my/new/scratch/dir"
