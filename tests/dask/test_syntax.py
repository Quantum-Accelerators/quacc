from __future__ import annotations

import pytest

dask = pytest.importorskip("dask")
pytest.importorskip("distributed")

from dask.distributed import get_client

from quacc import flow, job, subflow

client = get_client()


def test_dask_decorators(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @subflow
    def add_distributed2(vals, c, op):
        return [op(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def dynamic_workflow2(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed2(result2, c, add)

    @flow
    def dynamic_workflow3(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        result3 = add_distributed(result2, c)
        result4 = add_distributed(result3, c)
        return add(result4[0], c)

    assert client.compute(add(1, 2)).result() == 3
    assert client.compute(mult(1, 2)).result() == 2
    assert client.compute(workflow(1, 2, 3)).result() == 9
    assert client.compute(dynamic_workflow(1, 2, 3)).result() == [6, 6, 6]
    assert client.compute(dynamic_workflow2(1, 2, 3)).result() == [6, 6, 6]
    assert client.compute(dynamic_workflow3(1, 2, 3)).result() == 12


def test_dask_decorators_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @job()
    def make_more(val):
        return [val] * 3

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @subflow()
    def add_distributed2(vals, c, op):
        return [op(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow()
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow()
    def dynamic_workflow2(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed2(result2, c, add)

    assert client.compute(add(1, 2)).result() == 3
    assert client.compute(mult(1, 2)).result() == 2
    assert client.compute(workflow(1, 2, 3)).result() == 9
    assert client.compute(dynamic_workflow(1, 2, 3)).result() == [6, 6, 6]
    assert client.compute(dynamic_workflow2(1, 2, 3)).result() == [6, 6, 6]
