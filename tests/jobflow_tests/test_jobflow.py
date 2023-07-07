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
        if "quacc-tmp" in f or "quacc_" in f or "job_" in f or f == "tmp_dir":
            rmtree(f)


STORE = jf.JobStore(MemoryStore())


def test_tutorial1():
    from ase.build import bulk
    from jobflow import Flow, job, run_locally

    from quacc.recipes.emt.core import relax_job, static_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define Job 1
    job1 = job(relax_job)(atoms)

    # Define Job 2, which takes the output of Job 1 as input
    job2 = job(static_job)(job1.output)

    # Define the workflow
    workflow = Flow([job1, job2])

    # Run the workflow locally
    responses = run_locally(workflow, store=STORE, create_folders=True)

    # Get the result
    result = responses[job2.uuid][1].output


def test_tutorial2():
    from ase.build import bulk, molecule
    from jobflow import Flow, job, run_locally

    from quacc.recipes.emt.core import relax_job

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    job1 = job(relax_job)(atoms1)
    job2 = job(relax_job)(atoms2)

    # Define the workflow
    workflow = Flow([job1, job2])

    # Run the workflow locally
    responses = run_locally(workflow, store=STORE, create_folders=True)

    # Get the result
    result = responses[job2.uuid][1].output


def test_tutorial3():
    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(bulk_to_slabs_flow)(job1.output)
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
    flow = jf.Flow([job1, job2])

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

    store = jf.JobStore(MemoryStore())

    atoms = bulk("Cu")

    job = jf.job(bulk_to_slabs_flow)(atoms, slab_static_job=None)
    jf.run_locally(job, store=store, ensure_success=True, create_folders=True)

    job = jf.job(bulk_to_slabs_flow)(
        atoms,
        slab_static_job=None,
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_kwargs": {"asap_cutoff": True},
            "relax_cell": False,
        },
    )
    jf.run_locally(job, store=store, ensure_success=True, create_folders=True)

    job = jf.job(bulk_to_slabs_flow)(
        atoms,
        slabgen_kwargs={"max_slabs": 2},
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_kwargs": {"asap_cutoff": True},
            "relax_cell": False,
        },
    )
    responses = jf.run_locally(
        job, store=store, ensure_success=True, create_folders=True
    )

    assert len(responses) == 5
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    assert "generated_slabs" in output0
    assert len(output0["generated_slabs"][0]) == 64

    output1 = responses[uuids[1]][1].output
    assert output1["nsites"] == 64
    assert output1["parameters"]["asap_cutoff"] is True
    assert output1["name"] == "EMT Relax"

    output2 = responses[uuids[-1]][1].output
    assert output2["nsites"] == 80
    assert output2["parameters"]["asap_cutoff"] is False
    assert output2["name"] == "EMT Static"
