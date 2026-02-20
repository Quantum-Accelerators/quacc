from __future__ import annotations

import pytest

redun = pytest.importorskip("redun")
from ase.build import bulk, molecule

from quacc import flow
from quacc.recipes.emt.core import relax_job, static_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


@pytest.fixture
def scheduler():
    return redun.Scheduler()


def test_tutorial1a(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the job locally
    assert "atoms" in scheduler.run(relax_job(atoms))


def test_tutorial1b(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the job locally
    assert len(scheduler.run(bulk_to_slabs_flow(atoms))) == 8


def test_tutorial2a(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    @flow
    def workflow(atoms):
        result1 = relax_job(atoms)  # (1)!
        return static_job(result1["atoms"])  # (2)!

    atoms = bulk("Cu")

    # Define the workflow
    assert "atoms" in scheduler.run(workflow(atoms))


def test_tutorial2b(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    # Define workflow
    @flow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return {"result1": result1, "result2": result2}

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Dispatch the workflow
    assert "atoms" in scheduler.run(workflow(atoms1, atoms2))["result1"]


def test_tutorial2c(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk["atoms"], run_static=False)  # (1)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    assert len(scheduler.run(workflow(atoms))) == 4
