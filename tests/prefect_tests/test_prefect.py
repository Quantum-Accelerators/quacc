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


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    if prefect:
        with prefect_test_harness():
            yield


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial1(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk
    from prefect import flow, task
    from prefect.task_runners import SequentialTaskRunner

    from quacc.recipes.emt.core import relax_job, static_job

    # Define the workflow
    @flow(task_runner=SequentialTaskRunner())  # (3)!
    def workflow(atoms):
        # Call Task 1
        future1 = task(relax_job).submit(atoms)  # (4)!

        # Call Task 2, which takes the output of Task 1 as input
        future2 = task(static_job).submit(future1)

        return future2

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow with Prefect tracking
    workflow(atoms).result()


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial2(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk, molecule
    from prefect import flow, task
    from prefect.task_runners import SequentialTaskRunner

    from quacc.recipes.emt.core import relax_job

    # Define workflow
    @flow(task_runner=SequentialTaskRunner())
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        future1 = task(relax_job).submit(atoms1)
        future2 = task(relax_job).submit(atoms2)

        return future1, future2

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Run the workflow with Prefect tracking
    future1, future2 = workflow(atoms1, atoms2)
    future1.result(), future2.result()


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial3(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk
    from prefect import flow, task
    from prefect.task_runners import SequentialTaskRunner

    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    @flow(task_runner=SequentialTaskRunner())
    def workflow(atoms):
        future1 = task(relax_job).submit(atoms)
        future2 = task(bulk_to_slabs_flow).submit(future1, slab_static=None)

        return future2

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    workflow(atoms).result()


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial4(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk
    from prefect import flow, task
    from prefect.task_runners import SequentialTaskRunner

    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

    bulk_to_slabs_flow.task_runner = SequentialTaskRunner()

    @flow(task_runner=SequentialTaskRunner())
    def workflow(atoms):
        future1 = task(relax_job).submit(atoms)
        result = bulk_to_slabs_flow(future1, run_slab_static=False)

        return result

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    result = workflow(atoms)


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

    from prefect.task_runners import SequentialTaskRunner

    from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

    bulk_to_slabs_flow.task_runner = SequentialTaskRunner()
    bulk_to_slabs_flow(bulk("Cu"))
