from __future__ import annotations

import os
from pathlib import Path

from ase.build import bulk

from quacc import change_settings
from quacc.recipes.emt.core import relax_job
from quacc.settings import QuaccSettings

FILE_DIR = Path(__file__).parent


def test_file(tmp_path, monkeypatch):
    with open(tmp_path / "quacc_test.yaml", "w") as f:
        f.write("GZIP_FILES: false\nWORKFLOW_ENGINE: None")
    monkeypatch.setenv("QUACC_CONFIG_FILE", os.path.join(tmp_path, "quacc_test.yaml"))

    assert QuaccSettings().GZIP_FILES is False
    assert QuaccSettings().WORKFLOW_ENGINE is None
    os.remove(tmp_path / "quacc_test.yaml")


def test_results_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    output = relax_job(atoms)
    assert "opt.traj.gz" in os.listdir(output["dir_name"])

    with change_settings({"GZIP_FILES": False}):
        output = relax_job(atoms)
        assert "opt.traj" in os.listdir(output["dir_name"])


def test_env_var(monkeypatch, tmp_path):
    p = tmp_path / "my/scratch/dir"
    monkeypatch.setenv("QUACC_SCRATCH_DIR", p)
    assert p.expanduser().resolve() == QuaccSettings().SCRATCH_DIR


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

    monkeypatch.setenv("QUACC_GZIP_FILES", "False")
    assert QuaccSettings().GZIP_FILES is False

    monkeypatch.setenv("QUACC_GZIP_FILES", "false")
    assert QuaccSettings().GZIP_FILES is False


def test_yaml(tmp_path, monkeypatch):
    p = tmp_path / "my/new/scratch/dir"
    monkeypatch.delenv("QUACC_SCRATCH_DIR", raising=False)
    with open(tmp_path / "quacc_test.yaml", "w") as f:
        f.write(f"SCRATCH_DIR: {p}")
    monkeypatch.setenv("QUACC_CONFIG_FILE", os.path.join(tmp_path, "quacc_test.yaml"))
    assert p.expanduser().resolve() == QuaccSettings().SCRATCH_DIR
