import pytest
from ase.build import bulk

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
    assert "atoms" in scheduler.run(bulk_to_slabs_flow(atoms))[0]


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
