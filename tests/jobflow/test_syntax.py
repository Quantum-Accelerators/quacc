from __future__ import annotations

import warnings

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
    assert isinstance(workflow(1, 2, 3), jf.Flow)
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
    assert isinstance(workflow(1, 2, 3), jf.Flow)
    assert isinstance(add_distributed([1, 2, 3], 4), jf.Job)


def test_jobflow_tuple_argument(tmp_path, monkeypatch, jobflow_output):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def sum_values(values):
        return sum(values)

    @flow
    def workflow():
        job1 = add(1, 2)
        job2 = add(3, 4)
        return (sum_values((job1, job2)),)

    workflow_flow = workflow()
    responses = jf.run_locally(workflow_flow, ensure_success=True)
    assert jobflow_output(responses, workflow_flow.output[0]) == 10


def test_jobflow_flow_output_warning_is_suppressed(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @flow
    def workflow(a, b):
        return add(a, b)

    with warnings.catch_warnings(record=True) as caught_warnings:
        assert isinstance(workflow(1, 2), jf.Flow)

    assert not caught_warnings


def test_jobflow_dict_output(tmp_path, monkeypatch, jobflow_output):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @flow
    def workflow():
        return {"sum": add(1, 2)}

    workflow_flow = workflow()
    responses = jf.run_locally(workflow_flow, ensure_success=True)
    assert jobflow_output(responses, workflow_flow.output["sum"]) == 3
