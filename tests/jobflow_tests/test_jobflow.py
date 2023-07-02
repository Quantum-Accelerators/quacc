import os
from shutil import rmtree

import jobflow as jf
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job
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


def comparison1():
    @jf.job
    def add(a, b):
        return a + b

    @jf.job
    def mult(a, b):
        return a * b

    job1 = add(1, 2)
    job2 = mult(job1.output, 3)
    flow = jf.Flow([job1, job2], output=job2.output)

    responses = jf.run_locally(flow, ensure_success=True)
    assert responses[job2.uuid][1].output == 9


def comparison2():
    @jf.job
    def add(a, b):
        return a + b

    @jf.job
    def make_more(val):
        return [val] * 3

    @jf.job
    def add_distributed(vals, c):
        jobs = []
        for val in vals:
            jobs.append(add(val, c))

        flow = jf.Flow(jobs)
        return jf.Response(replace=flow)

    job1 = add(1, 2)
    job2 = make_more(job1.output)
    job3 = add_distributed(job2.output, 3)
    flow = jf.Flow([job1, job2, job3])

    jf.run_locally(flow, ensure_success=True)  # [6, 6, 6] in final 3 jobs


def test_emt_flow():
    from quacc.recipes.emt.jobflow.slabs import bulk_to_slabs_flow

    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(bulk_to_slabs_flow)(job1.output["atoms"])
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, store=STORE, create_folders=True, ensure_success=True)
