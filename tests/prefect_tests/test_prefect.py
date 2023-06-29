import os
from shutil import rmtree

import pytest
from ase.build import bulk, molecule

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

try:
    import prefect
    from prefect import flow, task
    from prefect.testing.utilities import prefect_test_harness

except ImportError:
    prefect = None


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or ".gz" in f
        ):
            os.remove(f)
        if "quacc-tmp" in f or "job_" in f or f == "tmp_dir" or f == "runinfo":
            rmtree(f)


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    if prefect:
        with prefect_test_harness():
            yield


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial1():
    # Define the workflow
    @flow
    def workflow(atoms):
        # Call Task 1
        result1 = task(relax_job)(atoms)

        # Call Task 2, which takes the output of Task 1 as input
        result2 = task(static_job)(result1["atoms"])

        return result2

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow with Prefect tracking
    workflow(atoms)


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial2():
    # Define workflow
    @flow
    def workflow2(atoms1, atoms2):
        # Define two independent relaxation jobs
        result1 = task(relax_job)(atoms1)
        result2 = task(relax_job)(atoms2)

        return {"result1": result1, "result2": result2}

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Run the workflow with Prefect tracking
    workflow2(atoms1, atoms2)


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial3():
    @flow
    def workflow(atoms):
        relaxed_bulk = task(relax_job)(atoms)
        relaxed_slabs = task(bulk_to_slabs_flow)(
            relaxed_bulk["atoms"], slab_static_electron=None
        )

        return relaxed_slabs

    atoms = bulk("Cu")
    workflow(atoms)


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial4():
    from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

    @flow
    def workflow(atoms):
        relaxed_bulk = task(relax_job)(atoms)
        relaxed_slabs = bulk_to_slabs_flow(
            relaxed_bulk["atoms"], slab_static_electron=None
        )

        return relaxed_slabs

    atoms = bulk("Cu")
    workflow(atoms)
