import pytest
from ase.build import bulk

from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()

SETTINGS.WORKFLOW_ENGINE = "prefect"

try:
    import prefect
    from prefect.testing.utilities import prefect_test_harness

except ImportError:
    prefect = None


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    if prefect:
        with prefect_test_harness():
            yield


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial1(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

    # Define the workflow
    @flow  # (1)!
    def workflow(atoms):
        # Call Task 1
        future1 = relax_job.submit(atoms)  # (2)!

        # Call Task 2, which takes the output of Task 1 as input
        future2 = static_job.submit(future1)

        return future2

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow with Prefect tracking
    result = workflow(atoms).result()  # (3)!
    assert "atoms" in result


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial2(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk, molecule

    from quacc import flow
    from quacc.recipes.emt.core import relax_job

    # Define workflow
    @flow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        future1 = relax_job.submit(atoms1)
        future2 = relax_job.submit(atoms2)

        return future1, future2

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Run the workflow with Prefect tracking
    future1, future2 = workflow(atoms1, atoms2)
    assert "atoms" in future1.result()
    assert "atoms" in future2.result()


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial3(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

    @flow
    def workflow(atoms):
        future1 = relax_job.submit(atoms)
        slab_futures = bulk_to_slabs_flow(future1, run_slab_static=False)  # (1)!

        return slab_futures

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    slab_futures = workflow(atoms)
    result = [slab_future.result() for slab_future in slab_futures]  # (2)!
    assert len(result) == 4


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_comparison1():
    from prefect import flow, task

    @task
    def add(a, b):
        return a + b

    @task
    def mult(a, b):
        return a * b

    @flow
    def workflow(a, b, c):
        return mult.submit(add.submit(a, b), c)

    assert workflow(1, 2, 3).result() == 9


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_comparison2():
    from prefect import flow, task

    @task
    def add(a, b):
        return a + b

    @flow
    def workflow(a, b, c):
        future1 = add.submit(a, b)

        vals_to_add = [future1.result()] * 3
        return [add.submit(val, c).result() for val in vals_to_add]

    assert workflow(1, 2, 3) == [6, 6, 6]


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_emt_flow(tmpdir):
    tmpdir.chdir()

    from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

    bulk_to_slabs_flow(bulk("Cu"))
