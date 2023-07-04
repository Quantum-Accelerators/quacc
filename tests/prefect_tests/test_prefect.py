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
        if "quacc-tmp" in f or "job_" in f or f == "tmp_dir":
            rmtree(f)


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    if prefect:
        with prefect_test_harness():
            yield


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial1():
    @flow
    def workflow(atoms):
        # Call Task 1
        future1 = task(relax_job).submit(atoms)

        # Call Task 2, which takes the output of Task 1 as input
        future2 = task(static_job).submit(future1.result())

        return future2.result()

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow with Prefect tracking
    workflow(atoms)


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial2():
    @flow
    def workflow2(atoms1, atoms2):
        # Define two independent relaxation jobs
        future1 = task(relax_job).submit(atoms1)
        future2 = task(relax_job).submit(atoms2)

        return {"result1": future1.result(), "result2": future2.result()}

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Run the workflow with Prefect tracking
    workflow2(atoms1, atoms2)


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial3():
    @flow
    def workflow(atoms):
        future1 = task(relax_job).submit(atoms)
        future2 = task(bulk_to_slabs_flow).submit(
            future1.result()["atoms"], slab_static_electron=None
        )

        return future2.result()

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    workflow(atoms)


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_comparison1():
    @task
    def add(a, b):
        return a + b

    @task
    def mult(a, b):
        return a * b

    @flow
    def workflow(a, b, c):
        future1 = add.submit(a, b)
        future2 = mult.submit(future1.result(), c)
        return future2

    assert workflow(1, 2, 3).result() == 9


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_comparison2():
    @task
    def add(a, b):
        return a + b

    @task
    def make_more(val):
        return [val] * 3

    @flow
    def workflow(a, b, c):
        future1 = add.submit(a, b)
        future2 = make_more.submit(future1.result())
        return [add.submit(val, c).result() for val in future2.result()]

    assert workflow(1, 2, 3) == [6, 6, 6]


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_emt_flow():
    from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

    @flow
    def workflow(atoms):
        future1 = task(relax_job).submit(atoms)
        result = bulk_to_slabs_flow(future1.result(), task(relax_job), task(static_job))

        return result

    atoms = bulk("Cu")
    workflow(atoms)

    @flow
    def workflow2(atoms):
        future1 = task(relax_job).submit(atoms)
        result = bulk_to_slabs_flow(future1.result()["atoms"], task(relax_job), None)

        return result

    atoms = bulk("Cu")
    workflow2(atoms)
