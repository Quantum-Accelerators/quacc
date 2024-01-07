import pytest

parsl = pytest.importorskip("parsl")

import contextlib

from ase.build import bulk, molecule

from quacc import job, subflow
from quacc.recipes.emt.core import relax_job, static_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


def setup_module():
    with contextlib.suppress(Exception):
        parsl.load()


def teardown_module():
    parsl.clear()


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call the PythonApp
    future = relax_job(atoms)  # (1)!

    # Print result
    assert "atoms" in future.result()  # (2)!


def test_tutorial1b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call the PythonApp
    future = relax_job(atoms)  # (1)!

    # Print result
    assert "atoms" in future.result()  # (2)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    future = bulk_to_slabs_flow(atoms)  # (1)!

    # Print the results
    assert "atoms" in future.result()[0]  # (2)!


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define the workflow
    def workflow(atoms):
        # Define Job 1
        future1 = relax_job(atoms)  # (1)!

        # Define Job 2, which takes the output of Job 1 as input
        return static_job(future1["atoms"])

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow
    future = workflow(atoms)

    # Fetch the result
    result = future.result()  # (2)!
    assert "atoms" in result


def test_tutorial2b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define workflow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return {"result1": result1, "result2": result2}

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    futures = workflow(atoms1, atoms2)

    # Fetch the results
    result1 = futures["result1"].result()
    result2 = futures["result2"].result()

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
    future = workflow(atoms)

    # Fetch the results
    result = future.result()

    # Print the results
    assert len(result) == 4


def test_comparison1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    def workflow(a, b, c):  #  (2)!
        return mult(add(a, b), c)

    result = workflow(1, 2, 3).result()  # 9
    assert result == 9


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

    future1 = add(1, 2)
    future2 = make_more(future1)
    future3 = add_distributed(future2, 3)

    assert future3.result() == [6, 6, 6]


def test_comparison3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    future1 = add(1, 2)
    future2 = mult(future1, 3)

    assert future2.result() == 9
