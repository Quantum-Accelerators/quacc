import os
from pathlib import Path

from ase.build import bulk
from maggma.stores import MemoryStore

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job
from quacc.settings import QuaccSettings

DEFAULT_SETTINGS = SETTINGS.model_copy()


FILE_DIR = Path(__file__).parent


def setup_function():
    SETTINGS.STORE = None
    SETTINGS.GZIP_FILES = True
    SETTINGS.CREATE_UNIQUE_DIR = False


def teardown_function():
    SETTINGS.STORE = DEFAULT_SETTINGS.STORE
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.CREATE_UNIQUE_DIR = DEFAULT_SETTINGS.CREATE_UNIQUE_DIR


def test_file_v1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with open("quacc_test.yaml", "w") as f:
        f.write("GZIP_FILES: false\nWORKFLOW_ENGINE:\nDEBUG: true")
    monkeypatch.setenv(
        "QUACC_CONFIG_FILE", os.path.join(os.getcwd(), "quacc_test.yaml")
    )

    assert QuaccSettings().GZIP_FILES is False
    assert QuaccSettings().WORKFLOW_ENGINE is None
    assert QuaccSettings().DEBUG is True
    os.remove("quacc_test.yaml")


def test_store(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    SETTINGS.STORE = MemoryStore()
    atoms = bulk("Cu")
    static_job(atoms)


def test_results_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    relax_job(atoms)
    assert "opt.traj.gz" in os.listdir(os.getcwd())
    os.remove("opt.traj.gz")
    SETTINGS.GZIP_FILES = False
    relax_job(atoms)
    assert "opt.traj" in os.listdir(os.getcwd())
    os.remove("opt.traj")


def test_env_var(monkeypatch, tmp_path):
    p = tmp_path / "my/scratch/dir"
    monkeypatch.setenv("QUACC_SCRATCH_DIR", p)
    assert p.expanduser().resolve() == QuaccSettings().SCRATCH_DIR


def test_yaml(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    p = tmp_path / "my/new/scratch/dir"
    monkeypatch.delenv("QUACC_SCRATCH_DIR", raising=False)
    with open("quacc_test.yaml", "w") as f:
        f.write(f"SCRATCH_DIR: {p}")
    monkeypatch.setenv("QUACC_CONFIG_FILE", "quacc_test.yaml")
    assert p.expanduser().resolve() == QuaccSettings().SCRATCH_DIR
