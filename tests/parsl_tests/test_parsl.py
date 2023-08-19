import pytest

from quacc import SETTINGS

try:
    import parsl
    from parsl import join_app, python_app
except ImportError:
    parsl = None

WFLOW_ENGINE = SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None


@pytest.mark.skipif(
    parsl is None or WFLOW_ENGINE != "parsl",
    reason="Parsl is not installed or specified in config",
)
def test_tutorial1(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job, static_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call App 1
    future1 = relax_job(atoms)  # (1)!

    # Call App 2, which takes the output of App 1 as input
    future2 = static_job(future1)

    assert "atoms" in future2.result()


@pytest.mark.skipif(
    parsl is None or WFLOW_ENGINE != "parsl",
    reason="Parsl is not installed or specified in config",
)
def test_tutorial2(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk, molecule

    from quacc.recipes.emt.core import relax_job

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    future1 = relax_job(atoms1)
    future2 = relax_job(atoms2)

    # Print the results
    assert "atoms" in future1.result()
    assert "atoms" in future2.result()


@pytest.mark.skipif(
    parsl is None or WFLOW_ENGINE != "parsl",
    reason="Parsl is not installed or specified in config",
)
def test_tutorial3(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    future1 = relax_job(atoms)
    future2 = bulk_to_slabs_flow(future1, slab_static=None)  # (1)!

    assert len(future2.result()) == 4
    assert future2.done()


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed",
)
def test_comparison1(tmpdir):
    tmpdir.chdir()

    @python_app
    def add(a, b):
        return a + b

    @python_app
    def mult(a, b):
        return a * b

    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert workflow(1, 2, 3).result() == 9


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed",
)
def test_comparison2(tmpdir):
    tmpdir.chdir()

    @python_app
    def add(a, b):
        return a + b

    @python_app
    def make_more(val):
        return [val] * 3

    @join_app
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    def workflow(a, b, c):
        future1 = add(a, b)
        future2 = make_more(future1)
        return add_distributed(future2, c)

    assert workflow(1, 2, 3).result() == [6, 6, 6]
