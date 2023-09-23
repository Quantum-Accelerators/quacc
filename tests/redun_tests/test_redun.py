import pytest

from quacc import SETTINGS

try:
    import redun

    redun = redun if SETTINGS.WORKFLOW_ENGINE == "redun" else None
except ImportError:
    redun = None

DEFAULT_SETTINGS = SETTINGS.copy()


@pytest.mark.skipif(redun is None, reason="Redun not installed")
def test_redun_decorators(tmpdir):
    tmpdir.chdir()
    from redun import Scheduler

    SETTINGS.WORKFLOW_ENGINE = "redun"
    scheduler = Scheduler()

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

    assert scheduler.run(add(1, 2)) == 3
    assert scheduler.run(mult(1, 2)) == 2
    assert scheduler.run(workflow(1, 2, 3)) == 9
    assert scheduler.run(dynamic_workflow(1, 2, 3)) == [6, 6, 6]


@pytest.mark.skipif(
    redun is None,
    reason="Redun is not installed or specified in config",
)
def test_tutorial1a(tmpdir):
    tmpdir.chdir()
    from ase.build import bulk
    from redun import Scheduler

    from quacc.recipes.emt.core import relax_job

    scheduler = Scheduler()

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the job locally
    assert "atoms" in scheduler.run(relax_job(atoms))


@pytest.mark.skipif(
    redun is None,
    reason="Redun is not installed or specified in config",
)
def test_tutorial1b(tmpdir):
    tmpdir.chdir()
    from ase.build import bulk
    from redun import Scheduler

    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    scheduler = Scheduler()

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the job locally
    assert len(scheduler.run(bulk_to_slabs_flow(atoms))) == 4


@pytest.mark.skipif(
    redun is None,
    reason="Redun is not installed or specified in config",
)
def test_tutorial2a(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk
    from redun import Scheduler

    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

    scheduler = Scheduler()

    @flow
    def workflow(atoms):
        result1 = relax_job(atoms)  # (1)!
        return static_job(result1)  # (2)!

    atoms = bulk("Cu")

    # Define the workflow
    assert "atoms" in scheduler.run(workflow(atoms))


@pytest.mark.skipif(
    redun is None,
    reason="Redun is not installed or specified in config",
)
def test_tutorial2b(tmpdir):
    tmpdir.chdir()
    from ase.build import bulk, molecule
    from redun import Scheduler

    from quacc import flow
    from quacc.recipes.emt.core import relax_job

    # Instantiate the scheduler
    scheduler = Scheduler()

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
    assert "atoms" in scheduler.run(workflow(atoms1, atoms2))["result1"]


@pytest.mark.skipif(
    redun is None,
    reason="Redun is not installed or specified in config",
)
def test_tutorial2c(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk
    from redun import Scheduler

    from quacc import flow
    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    scheduler = Scheduler()

    # Define the workflow
    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk, run_static=False)  # (1)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    assert len(scheduler.run(workflow(atoms))) == 4
