from __future__ import annotations

import pytest

prefect = pytest.importorskip("prefect")

import time

import psutil

from quacc import change_settings, flow, job, subflow

n_cpus = psutil.cpu_count() or 1


def test_patch():
    @job
    def task1(val: float) -> dict:
        return {"input": val, "result": val * 100}

    @job
    def task2(val: float) -> dict:
        return {"input": val, "result": val * 200}

    @flow
    def workflow(val):
        future1 = task1(val)
        return task2(future1["result"])

    assert workflow(1) == {"input": 100, "result": 20000}


def test_patch_local():
    with change_settings({"PREFECT_AUTO_SUBMIT": False}):

        @job
        def task1(val: float) -> dict:
            return {"input": val, "result": val * 100}

        @job
        def task2(val: float) -> dict:
            return {"input": val, "result": val * 200}

        @flow
        def workflow(val):
            future1 = task1(val)
            return task2(future1["result"])

        assert workflow(1) == {"input": 100, "result": 20000}


def test_prefect_decorators1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @flow
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2) == 3


def test_prefect_decorators2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert workflow(1, 2, 3) == 9


def test_prefect_decorators3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert dynamic_workflow(1, 2, 3) == [6, 6, 6]


def test_prefect_decorators4(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed2(vals, c, op):
        return [op(val, c) for val in vals]

    @flow
    def dynamic_workflow2(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed2(result2, c, add)

    assert dynamic_workflow2(1, 2, 3) == [6, 6, 6]


def test_prefect_decorators5(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def dynamic_workflow3(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        result3 = add_distributed(result2, c)
        result4 = add_distributed(result3, c)
        return add(result4[0], c)

    assert dynamic_workflow3(1, 2, 3) == 12


def test_prefect_decorators_local(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with change_settings({"PREFECT_AUTO_SUBMIT": False}):

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

        @flow
        def add_flow(a, b):
            return add(a, b)

        assert add_flow(1, 2) == 3
        assert workflow(1, 2, 3) == 9
        results = dynamic_workflow(1, 2, 3)
        assert results == [6, 6, 6]
        results = dynamic_workflow2(1, 2, 3)
        assert results == [6, 6, 6]
        assert dynamic_workflow3(1, 2, 3) == 12


def test_state_patch():
    @job
    def testjob(n):
        return {"result": n + 5}

    @flow
    def subflow(n):
        return testjob(n + 1)

    @flow
    def my_flow():
        # Try running a job
        job_result = testjob(1)  # return State

        # Try using the job as an input to a flow
        flow_result = subflow(job_result["result"])

        # Try using the outputs of a flow as the input to a job
        testjob(flow_result["result"])

        return 2

    assert my_flow() == 2


@pytest.mark.skipif(n_cpus < 2, reason="Need at least 2 CPUs to test concurrency")
def test_concurrency():
    @job
    def add(a, b):
        time.sleep(5)
        return {"val": a + b}

    @flow
    def myflow():
        a = add(1, 2)["val"]
        b = add(1, 2)["val"]
        return a, b

    t0 = time.time()
    myflow()
    t1 = time.time()
    dt = t1 - t0
    assert dt < 10


async def test_prefect_decorators_async(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    async def add(a, b):
        return a + b

    @flow
    async def add_flow(a, b):
        return add(a, b)

    assert (await add_flow(1, 2)) == 3


async def test_prefect_decorators3_async(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    async def add(a, b):
        return a + b

    @job
    async def make_more(val):
        return [val] * 3

    @subflow
    async def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    async def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return await add_distributed(result2, c)

    assert (await dynamic_workflow(1, 2, 3)) == [6, 6, 6]
