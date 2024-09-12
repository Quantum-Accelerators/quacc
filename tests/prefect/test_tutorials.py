from __future__ import annotations

import pytest

prefect = pytest.importorskip("prefect")
from ase.build import bulk, molecule

from quacc import flow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the workflow
    @flow
    def workflow(atoms):
        return relax_job(atoms)  # (1)!

    # Dispatch the workflow
    result = workflow(atoms)  # (2)!

    assert "atoms" in result


def test_tutorial1b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow
    results = bulk_to_slabs_flow(atoms)  # (1)!

    # Print the results
    for result in results:
        assert "atoms" in result


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    @flow
    def workflow(atoms):
        # Call Task 1
        future1 = relax_job(atoms)

        # Call Task 2, which takes the output of Task 1 as input
        return static_job(future1["atoms"])

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow with Prefect tracking
    result = workflow(atoms)
    assert "atoms" in result


def test_tutorial2b(tmp_path, monkeypatch):
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
    futures = workflow(atoms1, atoms2)

    # Fetch the results
    result1 = futures["result1"]
    result2 = futures["result2"]
    assert "atoms" in result1
    assert "atoms" in result2


def test_tutorial2c(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk["atoms"], run_static=False)  # (1)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow and retrieve result
    results = workflow(atoms)
    assert len(results) == 4
    for result in results:
        assert "atoms" in result
