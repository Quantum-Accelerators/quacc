import pytest

from quacc import SETTINGS, flow, job, subflow

try:
    import prefect
    from prefect.testing.utilities import prefect_test_harness

    prefect = prefect if SETTINGS.WORKFLOW_ENGINE == "prefect" else None

except ImportError:
    prefect = None


@pytest.mark.skipif(
    prefect is None,
    reason="This test requires Prefect",
)
@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


@pytest.mark.skipif(
    prefect is None,
    reason="This test requires Prefect",
)
def test_prefect_decorators(tmpdir):
    tmpdir.chdir()

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2).result() == 3
    assert workflow(1, 2, 3).result() == 9
    results = dynamic_workflow(1, 2, 3)
    assert [result.result() for result in results] == [6, 6, 6]


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial1a(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job

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
    assert "atoms" in result


# @pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
# def test_tutorial1b(tmpdir):
#     tmpdir.chdir()
#     from ase.build import bulk

#     from quacc.recipes.emt.slabs import bulk_to_slabs_flow

#     # Define the Atoms object
#     atoms = bulk("Cu")

#     # Dispatch the workflow
#     futures = bulk_to_slabs_flow(atoms)  # (1)!

#     # Print the results
#     results = [future.result() for future in futures]
#     for result in results:
#         assert "atom" in result


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial2a(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

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
    assert "atoms" in result


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_tutorial2b(tmpdir):
    tmpdir.chdir()
    from ase.build import bulk, molecule

    from quacc import flow
    from quacc.recipes.emt.core import relax_job

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
    assert "atoms" in result1
    assert "atoms" in result2


# @pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
# def test_tutorial2c(tmpdir):
#     tmpdir.chdir()

#     from ase.build import bulk

#     from quacc import flow
#     from quacc.recipes.emt.core import relax_job
#     from quacc.recipes.emt.slabs import bulk_to_slabs_flow

#     # Define the workflow
#     @flow
#     def workflow(atoms):
#         relaxed_bulk = relax_job(atoms)
#         relaxed_slabs = bulk_to_slabs_flow(relaxed_bulk, run_static=False)  # (1)!

#         return relaxed_slabs

#     # Define the Atoms object
#     atoms = bulk("Cu")

#     # Dispatch the workflow and retrieve result
#     futures = workflow(atoms)
#     results = [future.result() for future in futures]
#     for result in results:
#         assert "atoms" in result


@pytest.mark.skipif(prefect is None, reason="Prefect is not installed")
def test_comparison(tmpdir):
    tmpdir.chdir()
    from quacc import flow, job

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow  #  (2)!
    def workflow(a, b, c):
        return mult(add(a, b), c)

    future = workflow(1, 2, 3)  # (3)!
    result = future.result()
    assert result == 9
