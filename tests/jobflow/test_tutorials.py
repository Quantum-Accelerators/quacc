from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")

from ase.build import bulk, molecule

from quacc import flow, job
from quacc.recipes.emt.core import relax_job, static_job  # skipcq: PYL-C0412


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the Job
    job = relax_job(atoms)  # (1)!

    # Run the job locally
    jf.run_locally(job, create_folders=True, ensure_success=True)


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define Job 1
    job1 = relax_job(atoms)  # (1)!

    # Define Job 2, which takes the output of Job 1 as input
    job2 = static_job(job1["atoms"])  # (2)!

    # Define the workflow
    workflow = jf.Flow([job1, job2])

    # Run the workflow locally
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


def test_tutorial2a_flow_decorator(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the workflow
    @flow
    def workflow(atoms):
        # Define Job 1
        job1 = relax_job(atoms)
        # Define Job 2, which takes the output of Job 1 as input
        static_job(job1["atoms"])

    # Run the workflow locally
    jf.run_locally(workflow(atoms), create_folders=True, ensure_success=True)


def test_tutorial2b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

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


def test_tutorial2b_flow_decorator(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define the workflow
    @flow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        relax_job(atoms1)
        relax_job(atoms2)

    # Run the workflow locally
    jf.run_locally(workflow(atoms1, atoms2), create_folders=True, ensure_success=True)


def test_job_getitem():
    @job
    def greetings(s):
        return {"hello": f"Hello {s}", "bye": f"Goodbye {s}"}

    @job
    def upper(s):
        return s.upper()

    @flow
    def greet(s):
        job1 = greetings(s)
        job2 = upper(job1["hello"])  # No need for `job1.output["hello"]`
        return job2.output

    workflow = greet("World")
    response = jf.run_locally(workflow)
    assert response[workflow.output.uuid][1].output == "HELLO WORLD"


def test_comparison1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

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


def test_comparison1_flow_decorator(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job  # (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow
    def workflow():
        job1 = add(1, 2)
        job2 = mult(job1.output, 3)
        return job2.output  # or `return job`

    f = workflow()
    response = jf.run_locally(f, ensure_success=True)
    assert response[f.output.uuid][1].output == 9


def test_comparison2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @jf.job
    def add(a, b):
        return a + b

    @jf.job
    def make_more(val):
        return [val] * 3

    @jf.job
    def add_distributed(vals, c):
        jobs = []
        jobs = [add(val, c) for val in vals]

        flow = jf.Flow(jobs)
        return jf.Response(replace=flow)

    job1 = add(1, 2)
    job2 = make_more(job1.output)
    job3 = add_distributed(job2.output, 3)
    flow = jf.Flow([job1, job2, job3])

    jf.run_locally(flow, ensure_success=True)  # [6, 6, 6] in final 3 jobs


def test_comparison2_flow_decorator(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @jf.job
    def add(a, b):
        return a + b

    @jf.job
    def make_more(val):
        return [val] * 3

    @jf.job
    def add_distributed(vals, c):
        jobs = [add(val, c) for val in vals]

        flow = jf.Flow(jobs)
        return jf.Response(replace=flow)

    @jf.flow
    def workflow():
        job1 = add(1, 2)
        job2 = make_more(job1.output)
        add_distributed(job2.output, 3)

    jf.run_locally(workflow(), ensure_success=True)  # [6, 6, 6] in final 3 jobs


def test_comparison3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

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


def test_comparison3_flow_decorator(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow
    def workflow():
        job1 = add(1, 2)
        mult(job1.output, 3)

    jf.run_locally(workflow(), ensure_success=True)
