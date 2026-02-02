from __future__ import annotations

import os

import pytest

jobflow = pytest.importorskip("jobflow")

from ase.build import bulk

from quacc import flow, job
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


@pytest.mark.parametrize("job_decorators", [None, {"relax_job": job()}])
def test_functools(tmp_path, monkeypatch, job_decorators):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    flow = bulk_to_slabs_flow(
        atoms,
        run_static=False,
        job_params={"relax_job": {"opt_params": {"fmax": 0.1}}},
        job_decorators=job_decorators,
    )
    jobflow.run_locally(flow, ensure_success=True)


def test_copy_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    @flow
    def myflow(atoms):
        result1 = relax_job(atoms)
        return relax_job(result1["atoms"], copy_files={result1["dir_name"]: "opt.*"})

    output = jobflow.run_locally(myflow(atoms))
    first_output = next(iter(output.values()))[1].output
    assert "atoms" in first_output


def test_folders(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    job = relax_job(atoms)
    jobflow.run_locally(job, ensure_success=True, create_folders=True)
    files = os.listdir(tmp_path)
    assert len(files) == 1
    assert files[0].startswith("job")
    assert "opt.log.gz" in os.listdir(tmp_path / files[0])


def test_relax_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    @flow
    def relax_flow(atoms):
        result1 = relax_job(atoms)
        return relax_job(result1["atoms"])

    jobflow.run_locally(relax_flow(atoms), ensure_success=True)


def test_relaxed_slabs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk["atoms"], run_static=False)

    jobflow.run_locally(workflow(atoms), ensure_success=True)
