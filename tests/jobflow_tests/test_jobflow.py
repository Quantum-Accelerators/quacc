import pytest
from maggma.stores import MemoryStore

from quacc import SETTINGS, flow, job, subflow

try:
    import jobflow as jf

    jf = jf if SETTINGS.WORKFLOW_ENGINE == "jobflow" else None
except ImportError:
    jf = None

if jf:
    STORE = jf.JobStore(MemoryStore())


@pytest.mark.skipif(jf is None, reason="Jobflow not installed")
def test_jobflow_decorators(tmpdir):
    tmpdir.chdir()

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert not isinstance(add, jf.Job)
    assert not isinstance(mult, jf.Job)
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
    assert isinstance(workflow(1, 2, 3), jf.Job)
    assert isinstance(add_distributed([1, 2, 3], 4)[0], jf.Job)


@pytest.mark.skipif(jf is None, reason="Jobflow not installed")
def test_jobflow_decorators_args(tmpdir):
    tmpdir.chdir()

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert not isinstance(add, jf.Job)
    assert not isinstance(mult, jf.Job)
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
    assert isinstance(workflow(1, 2, 3), jf.Job)
    assert isinstance(add_distributed([1, 2, 3], 4)[0], jf.Job)


@pytest.mark.skipif(
    jf is None,
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
    jf is None,
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
    jf is None,
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
    jf is None,
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
    jf is None,
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
    jf is None,
    reason="Jobflow is not installed or specified in config",
)
def test_comparison3(tmpdir):
    tmpdir.chdir()
    import jobflow as jf

    from quacc import job

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    job1 = add(1, 2)
    job2 = mult(job1.output, 3)
    flow = jf.Flow([job1, job2])

    jf.run_locally(flow, ensure_success=True)


@pytest.mark.skipif(
    jf is None,
    reason="Jobflow is not installed or specified in config",
)
def test_comparison4(tmpdir):
    tmpdir.chdir()

    import jobflow as jf

    from quacc import job

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @job
    def add_distributed(vals, c):
        jobs = []
        for val in vals:
            jobs.append(add(val, c))
        return jf.Response(replace=jf.Flow(jobs))

    job1 = add(1, 2)
    job2 = make_more(job1.output)
    job3 = add_distributed(job2.output, 3)
    flow = jf.Flow([job1, job2, job3])

    jf.run_locally(flow, ensure_success=True)
