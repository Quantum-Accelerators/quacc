from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")

from quacc import flow, job, subflow


def test_jobflow_decorators(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

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
    assert hasattr(add, "original")
    assert hasattr(mult, "original")
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
    workflow_flow = workflow(1, 2, 3)
    assert isinstance(workflow_flow, jf.Flow)
    responses = jf.run_locally(workflow_flow, ensure_success=True)
    assert responses[workflow_flow.jobs[-1].uuid][1].output == 9
    assert isinstance(add_distributed([1, 2, 3], 4), jf.Job)


def test_jobflow_decorators_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

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
    assert hasattr(add, "original")
    assert hasattr(mult, "original")
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
    workflow_flow = workflow(1, 2, 3)
    assert isinstance(workflow_flow, jf.Flow)
    responses = jf.run_locally(workflow_flow, ensure_success=True)
    assert responses[workflow_flow.jobs[-1].uuid][1].output == 9
    assert isinstance(add_distributed([1, 2, 3], 4), jf.Job)


def test_jobflow_nested_job_arguments(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def combine(values, factor):
        return sum(values["items"]) * factor

    @flow
    def inner_workflow():
        return add(5, 6)

    @flow
    def workflow():
        values = {"items": [add(1, 2), add(3, 4), inner_workflow()]}
        return combine(values, factor=add(2, 4))

    workflow_flow = workflow()
    responses = jf.run_locally(workflow_flow, ensure_success=True)
    assert responses[workflow_flow.jobs[-1].uuid][1].output == 126
