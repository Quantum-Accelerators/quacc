import pytest

from quacc import SETTINGS, flow, job, strip_decorator, subflow

dask = pytest.importorskip("dask")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "dask",
    reason="This test requires the Dask workflow engine",
)
from dask.distributed import default_client

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
    def dynamic_workflow2(a, b, c, op):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed2(result2, c, op)

    assert client.compute(add(1, 2)).result() == 3
    assert client.compute(mult(1, 2)).result() == 2
    assert client.compute(workflow(1, 2, 3)).result() == 9
    assert client.gather(client.compute(dynamic_workflow(1, 2, 3))) == [6, 6, 6]
    assert client.gather(client.compute(dynamic_workflow2(1, 2, 3, add))) == [6, 6, 6]


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

    @subflow
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
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def dynamic_workflow2(a, b, c, op):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed2(result2, c, op)

    assert client.compute(add(1, 2)).result() == 3
    assert client.compute(mult(1, 2)).result() == 2
    assert client.compute(workflow(1, 2, 3)).result() == 9
    assert client.gather(client.compute(dynamic_workflow(1, 2, 3))) == [6, 6, 6]
    assert client.gather(client.compute(dynamic_workflow2(1, 2, 3, add))) == [6, 6, 6]


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
