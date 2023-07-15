import os
from shutil import rmtree

import covalent as ct
import pytest
from ase.build import bulk, molecule

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


def teardown_function():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or ".gz" in f
        ):
            os.remove(f)
        if "quacc-tmp" in f or f == "tmp_dir":
            rmtree(f)


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_tutorials():
    # Test of the various tutorials

    # Quick start -------------------------------------------------
    workflow_start = ct.lattice(relax_job)
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow_start)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    @ct.lattice(executor="local")
    def workflow_start2(atoms):
        relaxed_bulk = relax_job(atoms)
        relaxed_slabs = bulk_to_slabs_flow(relaxed_bulk)
        return relaxed_slabs

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow_start2)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # Tutorials ---------------------------------------------------

    # Define the workflow
    @ct.lattice
    def workflow(atoms):
        return relax_job(atoms)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow to the Covalent server
    # with the bulk Cu Atoms object as the input
    dispatch_id = ct.dispatch(workflow)(atoms)

    # Fetch the result from the server
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------

    # Define the workflow
    workflow = ct.lattice(relax_job)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow to the Covalent server
    # with the bulk Cu Atoms object as the input
    dispatch_id = ct.dispatch(workflow)(atoms)

    # Fetch the result from the server
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------

    @ct.lattice
    def workflow1(atoms):
        result1 = relax_job(atoms)
        result2 = static_job(result1)
        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow1)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow2(atoms1, atoms2):
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return [result1, result2]

    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")
    dispatch_id = ct.dispatch(workflow2)(atoms1, atoms2)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow3(atoms):
        relaxed_bulk = relax_job(atoms)
        relaxed_slabs = bulk_to_slabs_flow(relaxed_bulk, slab_static=None)
        return relaxed_slabs

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow3)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice(executor="local")
    def workflow4(atoms):
        result1 = relax_job(atoms)
        result2 = static_job(result1)
        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow4)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.electron
    def relax_electron(atoms):
        return relax_job(atoms)

    @ct.electron
    def static_electron(atoms):
        return static_job(atoms)

    @ct.lattice
    def workflow5(atoms):
        relax_electron.executor = "dask"
        static_electron.executor = "local"
        output1 = relax_electron(atoms)
        output2 = static_electron(output1)
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow5)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_comparison1():
    @ct.electron
    def add(a, b):
        return a + b

    @ct.electron
    def mult(a, b):
        return a * b

    @ct.lattice
    def workflow(a, b, c):
        return mult(add(a, b), c)

    # Locally
    assert workflow(1, 2, 3) == 9

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_comparison2():
    @ct.electron
    def add(a, b):
        return a + b

    @ct.electron
    def make_more(val):
        return [val] * 3

    @ct.electron
    @ct.lattice
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @ct.lattice
    def workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    # Locally
    assert workflow(1, 2, 3) == [6, 6, 6]

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"
