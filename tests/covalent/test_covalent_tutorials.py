import pytest

ct = pytest.importorskip("covalent")

from ase.build import bulk, molecule

from quacc import flow, job, subflow
from quacc.recipes.emt.core import relax_job, static_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


def test_quickstart(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow
    dispatch_id = ct.dispatch(bulk_to_slabs_flow)(atoms)

    # Fetch the results
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the workflow
    workflow = flow(relax_job)  # (1)!

    # Dispatch the workflow to the Covalent server
    # with the bulk Cu Atoms object as the input
    dispatch_id = ct.dispatch(workflow)(atoms)  # (3)!

    # Fetch the result from the server
    result = ct.get_result(dispatch_id, wait=True)  # (4)!
    assert result.status == "COMPLETED"


def test_tutorial1b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(bulk_to_slabs_flow)(atoms)  # (1)!
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    @flow  # (1)!
    def workflow(atoms):
        # Define Job 1
        result1 = relax_job(atoms)  # (2)!

        # Define Job 2, which takes the output of Job 1 as input
        return static_job(result1["atoms"])

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow to the Covalent server
    # with the bulk Cu Atoms object as the input
    dispatch_id = ct.dispatch(workflow)(atoms)  # (3)!

    # Fetch the result from the server
    result = ct.get_result(dispatch_id, wait=True)  # (4)!
    assert result.status == "COMPLETED"


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

    # Dispatch the workflow to the Covalent server
    dispatch_id = ct.dispatch(workflow)(atoms1, atoms2)

    # Fetch the results from the server
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


def test_tutorial2c(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk["atoms"], run_static=False)  # (1)!

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


def test_tutorial_excecutor1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @flow(executor="local")
    def workflow4(atoms):
        result1 = relax_job(atoms)
        return static_job(result1["atoms"])

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow4)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


def test_tutorial_excecutor2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    relax_job.electron_object.executor = "dask"
    static_job.electron_object.executor = "local"

    @flow
    def workflow5(atoms):
        output1 = relax_job(atoms)
        return static_job(output1["atoms"])

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow5)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


def test_comparison1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job  # (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow  # (2)!
    def workflow(a, b, c):
        return mult(add(a, b), c)

    # Locally
    result = workflow(1, 2, 3)  # 9  (3)!

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)  # (4)!
    result = ct.get_result(dispatch_id, wait=True)  # 9  (5)!
    assert result.status == "COMPLETED"


def test_comparison2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow  # (1)!
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    # Locally
    result = workflow(1, 2, 3)  # e.g. [6, 6, 6]

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)
    result = ct.get_result(dispatch_id, wait=True)  # e.g. [6, 6, 6]

    assert result.status == "COMPLETED"


def test_comparison3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow  #  (2)!
    def workflow(a, b, c):
        return mult(add(a, b), c)

    dispatch_id = ct.dispatch(workflow)(1, 2, 3)  # (3)!
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"
