from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")
import os

from ase.build import bulk

from quacc.recipes.emt.core import relax_job


def test_folders(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    job = relax_job(atoms)
    jf.run_locally(job, ensure_success=True, create_folders=True)
    files = os.listdir(tmp_path)
    assert len(files) == 1
    assert "opt.log.gz" in os.listdir(tmp_path / files[0])
