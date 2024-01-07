import pytest

dask = pytest.importorskip("dask")

from dask.distributed import default_client

from quacc import flow, job, strip_decorator, subflow
from quacc.wflow_tools.customizers import update_parameters

client = default_client()


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

    assert client.compute(add(1, 2)).result() == 3
    assert client.compute(mult(1, 2)).result() == 2
    assert client.compute(workflow(1, 2, 3)).result() == 9
    assert client.compute(dynamic_workflow(1, 2, 3)).result() == [6, 6, 6]
    assert client.compute(dynamic_workflow2(1, 2, 3)).result() == [6, 6, 6]


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


def test_strip_decorators():
    @job
    def add(a, b):
        return a + b

    @flow
    def add2(a, b):
        return a + b

    @subflow
    def add3(a, b):
        return a + b

    stripped_add = strip_decorator(add)
    assert stripped_add(1, 2) == 3

    stripped_add2 = strip_decorator(add2)
    assert stripped_add2(1, 2) == 3

    stripped_add3 = strip_decorator(add3)
    assert stripped_add3(1, 2) == 3


def test_customize_funcs(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b=1):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c, d=0):
        return [add(val, c) + d for val in vals]

    @flow
    def dynamic_workflow(a, b, c=1):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def test_dynamic_workflow(a, b, c=3):
        result1 = add(a, b)
        result2 = make_more(result1)
        return update_parameters(add_distributed, {"d": 1}, decorator="subflow")(
            result2, c
        )

    add_ = update_parameters(add, {"b": 3}, decorator="job")
    dynamic_workflow_ = update_parameters(dynamic_workflow, {"c": 4}, decorator="flow")
    assert client.compute(add_(1)).result() == 4
    assert client.compute(dynamic_workflow_(1, 2)).result() == [7, 7, 7]
    assert client.compute(test_dynamic_workflow(1, 2)).result() == [7, 7, 7]
