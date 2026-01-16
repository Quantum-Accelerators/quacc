from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")
import os

from ase.build import bulk

from quacc.recipes.emt.core import relax_job
from quacc.wflow_tools.job_argument import Copy


def test_folders(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    job = relax_job(atoms)
    jf.run_locally(job, ensure_success=True, create_folders=True)
    files = os.listdir(tmp_path)
    assert len(files) == 1
    assert files[0].startswith("job")
    assert "opt.log.gz" in os.listdir(tmp_path / files[0])


def test_copy_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    job1 = relax_job(atoms)
    job2 = relax_job(
        job1.output["atoms"],
        copy_files=Copy(src_dir=job1.output["dir_name"], files="opt.*"),
    )
    flow = jf.Flow([job1, job2])
    jf.run_locally(flow, ensure_success=True, create_folders=True)
