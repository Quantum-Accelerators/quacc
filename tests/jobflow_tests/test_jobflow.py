import pytest
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc import SETTINGS

try:
    import jobflow as jf
except ImportError:
    jf = None

if jf:
    STORE = jf.JobStore(MemoryStore())
WFLOW_ENGINE = SETTINGS.WORKFLOW_ENGINE


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_tutorial1a(tmpdir):
    tmpdir.chdir()
    import jobflow as jf
    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the Job
    job = relax_job(atoms)  # (1)!

    # Run the job locally
    jf.run_locally(job, create_folders=True, ensure_success=True)


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_tutorial1b(tmpdir):
    tmpdir.chdir()
    import jobflow as jf
    from ase.build import bulk

    from quacc.recipes.emt._jobflow.slabs import bulk_to_slabs_flow

    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    flow = bulk_to_slabs_flow(atoms)

    # Run the workflow locally
    jf.run_locally(flow, create_folders=True, ensure_success=True)


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_tutorial2a(tmpdir):
    tmpdir.chdir()

    import jobflow as jf
    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job, static_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define Job 1
    job1 = relax_job(atoms)  # (1)!

    # Define Job 2, which takes the output of Job 1 as input
    job2 = static_job(job1.output)  # (2)!

    # Define the workflow
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_tutorial2b(tmpdir):
    tmpdir.chdir()

    import jobflow as jf
    from ase.build import bulk, molecule

    from quacc.recipes.emt.core import relax_job

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    job1 = relax_job(atoms1)
    job2 = relax_job(atoms2)

    # Define the workflow
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_tutorial2c(tmpdir):
    tmpdir.chdir()

    import jobflow as jf
    from ase.build import bulk

    from quacc.recipes.emt._jobflow.slabs import bulk_to_slabs_flow
    from quacc.recipes.emt.core import relax_job

    # Define the Atoms object
    atoms = bulk("Cu")

    # Construct the Flow
    job1 = relax_job(atoms)
    job2 = bulk_to_slabs_flow(job1.output, slab_static=None)
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_comparison1(tmpdir):
    tmpdir.chdir()

    import jobflow as jf

    from quacc import job

    @job  # (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    job1 = add(1, 2)
    job2 = mult(job1.output, 3)
    flow = jf.Flow([job1, job2])  # (2)!

    responses = jf.run_locally(flow, ensure_success=True)
    assert responses[job2.uuid][1].output == 9


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_comparison2(tmpdir):
    tmpdir.chdir()

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


@pytest.mark.skipif(
    jf is None or WFLOW_ENGINE != "jobflow",
    reason="Jobflow is not installed or specified in config",
)
def test_emt_flow(tmpdir):
    tmpdir.chdir()

    from quacc.recipes.emt._jobflow.slabs import bulk_to_slabs_flow

    store = jf.JobStore(MemoryStore())

    atoms = bulk("Cu")

    job = bulk_to_slabs_flow(
        atoms,
        slab_static=None,
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_swaps": {"asap_cutoff": True},
            "relax_cell": False,
        },
    )
    jf.run_locally(job, store=store, ensure_success=True, create_folders=True)

    job = bulk_to_slabs_flow(
        atoms,
        make_slabs_kwargs={"max_slabs": 2},
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_swaps": {"asap_cutoff": True},
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
