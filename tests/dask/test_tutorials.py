from __future__ import annotations

import pytest

dask = pytest.importorskip("dask")
pytest.importorskip("distributed")

from ase.build import bulk, molecule
from dask.distributed import get_client

from quacc.recipes.emt.core import relax_job, static_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412

client = get_client()


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call the PythonApp
    delayed = relax_job(atoms)  # (1)!

    # Print result
    assert "atoms" in client.compute(delayed).result()  # (2)!
    assert "atoms" in dask.compute(delayed)[0]


def test_tutorial1b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    delayed = bulk_to_slabs_flow(atoms)  # (1)!

    # Print the results
    assert "atoms" in client.compute(delayed).result()[0]


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    def workflow(atoms):
        # Define Job 1
        delayed1 = relax_job(atoms)  # (1)!

        # Define Job 2, which takes the output of Job 1 as input
        return static_job(delayed1["atoms"])

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow
    delayed = workflow(atoms)

    # Fetch the result
    assert "atoms" in client.compute(delayed).result()  # (2)!


def test_tutorial2b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define workflow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return [result1, result2]

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    delayed = workflow(atoms1, atoms2)

    # Fetch the results
    results = client.gather(client.compute(delayed))
    result1 = results[0]
    result2 = results[1]

    # Print the results
    assert "atoms" in result1
    assert "atoms" in result2


def test_tutorial2c(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk["atoms"], run_static=False)  # (1)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow
    delayed = workflow(atoms)

    # Fetch the results
    result = client.compute(delayed).result()

    # Print the results
    assert len(result) == 4
    assert "atoms" in result[0]
