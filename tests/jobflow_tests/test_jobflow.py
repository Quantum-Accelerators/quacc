import os
from shutil import rmtree

import jobflow as jf
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.jobflow.slabs import BulkToSlabsFlow
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


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


STORE = jf.JobStore(MemoryStore())


def test_tutorial1():
    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the compute job
    job = jf.job(static_job)(atoms)

    # Run the job locally
    jf.run_locally(job, store=STORE, create_folders=True, ensure_success=True)


def test_tutorial2():
    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the compute jobs
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    # Define the workflow
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, store=STORE, create_folders=True)


def test_tutorial3():
    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(bulk_to_slabs_flow)(job1.output["atoms"])
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, store=STORE, create_folders=True)


def test_emt_flow():
    atoms = bulk("Cu")

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
    jf.run_locally(workflow, store=STORE, create_folders=True, ensure_success=True)
