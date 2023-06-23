import os
from shutil import rmtree

import pytest
from ase.build import bulk, molecule
from parsl import python_app

from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

try:
    import parsl
except ImportError:
    parsl = None


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def setup_module():
    parsl.load()


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
        if "quacc-tmp-" in f or "job_" in f or f == "tmp_dir" or f == "runinfo":
            rmtree(f)


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial1():
    # Define the Python apps
    @python_app
    def relax_app(atoms):
        # All dependencies must be inside the Python app
        from quacc.recipes.emt.core import relax_job

        return relax_job(atoms)

    @python_app
    def static_app(atoms):
        # All dependencies must be inside the Python app
        from quacc.recipes.emt.core import static_job

        return static_job(atoms)

    # Define the workflow
    def workflow(atoms):
        # Call Job 1
        future1 = relax_app(atoms)

        # Call Job 2, which takes the output of Job 1 as input
        future2 = static_app(future1.result()["atoms"])

        return future2

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Run the workflow
    wf_future = workflow(atoms)
    wf_future.result()
    assert wf_future.done()


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial2():
    # Define the Python app
    @python_app
    def relax_app(atoms):
        # All dependencies must be inside the Python app
        from quacc.recipes.emt.core import relax_job

        return relax_job(atoms)

    # Define workflow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        future1 = relax_app(atoms1)
        future2 = relax_app(atoms2)

        return future1, future2

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Run the workflow
    future1, future2 = workflow(atoms1, atoms2)
    future1.result(), future2.result()
    assert future1.done()
    assert future2.done()


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial3():
    from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

    wf_future = bulk_to_slabs_flow(bulk("Cu"), slab_static_app=None)
    wf_future.result()
    assert wf_future.done()
