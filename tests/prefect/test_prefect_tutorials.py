prefect = pytest.importorskip("prefect")
import pytest
from ase.build import bulk, molecule
from prefect.testing.utilities import prefect_test_harness

from quacc import flow, job
from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the workflow
    @flow
    def workflow(atoms):
        return relax_job(atoms)  # (1)!

    # Dispatch the workflow
    future = workflow(atoms)  # (2)!

    # Fetch the result
    result = future.result()  # (3)!
    assert result.is_completed()


def test_tutorial1b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow
    futures = bulk_to_slabs_flow(atoms)  # (1)!

    # Print the results
    results = [future.result() for future in futures]
    for result in results:
        assert result.is_completed()


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    @flow
    def workflow(atoms):
        # Call Task 1
        future1 = relax_job(atoms)

        # Call Task 2, which takes the output of Task 1 as input
        return static_job(future1)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow with Prefect tracking
    result = workflow(atoms).result()
    assert result.is_completed()


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
    result1 = futures["result1"].result()
    result2 = futures["result2"].result()
    assert result1.is_completed()
    assert result2.is_completed()


def test_tutorial2c(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk, run_static=False)  # (1)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow and retrieve result
    futures = workflow(atoms)
    results = [future.result() for future in futures]
    for result in results:
        assert result.is_completed()
