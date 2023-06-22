import os
from shutil import rmtree

import jobflow as jf
import pytest
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job

try:
    import fireworks
except ImportError:
    fireworks = None


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


def test_emt():
    store = jf.JobStore(MemoryStore())

    atoms = bulk("Cu")

    # Test individual jobs
    job = jf.job(static_job)(atoms)
    jf.run_locally(job, store=store, ensure_success=True)

    job = jf.job(relax_job)(atoms)
    jf.run_locally(job, store=store, ensure_success=True)

    # Test a Flow
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


def test_emt_flow():
    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    @jf.job
    def relax_func(atoms):
        return relax_job(atoms)

    @jf.job
    def bulk_to_slabs_func(atoms):
        return bulk_to_slabs_flow(atoms, slab_static_electron=None)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    job1 = relax_func(atoms)
    job2 = bulk_to_slabs_func(job1.output["atoms"])
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


def test_emt_flow2():
    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.jobflow.slabs import BulkToSlabsFlow

    @jf.job
    def relax_func(atoms):
        return relax_job(atoms)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    job1 = relax_func(atoms)
    job2 = BulkToSlabsFlow().make(job1.output["atoms"])
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


@pytest.mark.skipif(fireworks is None, reason="This test requires fireworks")
def test_fireworks():
    from jobflow.managers.fireworks import flow_to_workflow, job_to_firework

    atoms = bulk("Cu")

    # Test fireworks creation
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    job_to_firework(job1)
    flow_to_workflow(workflow)
