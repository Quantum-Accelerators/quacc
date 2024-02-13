import os
from pathlib import Path

from ase.build import bulk
from maggma.stores import MemoryStore

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job, static_job
from quacc.settings import QuaccSettings

FILE_DIR = Path(__file__).parent


def test_file(tmp_path, monkeypatch):
    with open(tmp_path / "quacc_test.yaml", "w") as f:
        f.write("GZIP_FILES: false\nWORKFLOW_ENGINE: None\nDEBUG: True\nSTORE: null")
    monkeypatch.setenv("QUACC_CONFIG_FILE", os.path.join(tmp_path, "quacc_test.yaml"))

    assert QuaccSettings().GZIP_FILES is False
    assert QuaccSettings().WORKFLOW_ENGINE is None
    assert QuaccSettings().DEBUG is True
    assert QuaccSettings().STORE is None
    os.remove(tmp_path / "quacc_test.yaml")


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
    assert QuaccSettings().SCRATCH_DIR == p.expanduser().resolve()


def test_env_var2(monkeypatch, tmp_path):
    with open(tmp_path / "quacc_test.yaml", "w") as f:
        f.write("")

    monkeypatch.setenv("QUACC_CONFIG_FILE", os.path.join(tmp_path, "quacc_test.yaml"))
    assert str(QuaccSettings().CONFIG_FILE) == str(
        Path(os.path.join(tmp_path, "quacc_test.yaml"))
    )

    monkeypatch.setenv("QUACC_WORKFLOW_ENGINE", "None")
    assert QuaccSettings().WORKFLOW_ENGINE is None

    monkeypatch.setenv("QUACC_WORKFLOW_ENGINE", "null")
    assert QuaccSettings().WORKFLOW_ENGINE is None

    monkeypatch.setenv("QUACC_WORKFLOW_ENGINE", "none")
    assert QuaccSettings().WORKFLOW_ENGINE is None

    monkeypatch.setenv("QUACC_GZIP_FILES", "FaLsE")
    assert QuaccSettings().GZIP_FILES is False


def test_yaml(tmp_path, monkeypatch):
    p = tmp_path / "my/new/scratch/dir"
    monkeypatch.delenv("QUACC_SCRATCH_DIR", raising=False)
    with open(tmp_path / "quacc_test.yaml", "w") as f:
        f.write(f"SCRATCH_DIR: {p}")
    monkeypatch.setenv("QUACC_CONFIG_FILE", os.path.join(tmp_path, "quacc_test.yaml"))
    assert QuaccSettings().SCRATCH_DIR == p.expanduser().resolve()
