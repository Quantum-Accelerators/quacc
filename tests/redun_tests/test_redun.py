import pytest

from quacc import SETTINGS

try:
    import redun
except ImportError:
    redun = None

DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    SETTINGS.WORKFLOW_ENGINE = "redun"


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


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
        relaxed_slabs = bulk_to_slabs_flow(relaxed_bulk, slab_static=None)  # (1)!

        return relaxed_slabs

    # Define the Atoms object
    atoms = bulk("Cu")

    # Run the workflow
    assert len(scheduler.run(workflow(atoms))) == 4
