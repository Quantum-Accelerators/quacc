import contextlib

import pytest

from quacc import SETTINGS, flow, job, subflow

parsl = pytest.importorskip("dask")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "dask",
    reason="This test requires the Dask workflow engine",
)


def setup_module():
    from dask.distributed import Client

    Client()


def teardown_module():
    from dask.distributed import default_client

    default_client().close()


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

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2).compute() == 3
    assert mult(1, 2).compute() == 2
    assert workflow(1, 2, 3).compute() == 9
    assert dynamic_workflow(1, 2, 3).compute() == [6, 6, 6]


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

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow()
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2).compute() == 3
    assert mult(1, 2).compute() == 2
    assert workflow(1, 2, 3).compute() == 9
    assert dynamic_workflow(1, 2, 3).compute() == [6, 6, 6]
