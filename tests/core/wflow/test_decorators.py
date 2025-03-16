from __future__ import annotations

from maggma.stores import MemoryStore

from quacc import flow, job, subflow


def test_decorators(tmp_path, monkeypatch):
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

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert workflow(1, 2, 3) == 9
    assert dynamic_workflow(1, 2, 3) == [6, 6, 6]


def test_decorators_v2(tmp_path, monkeypatch):
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

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert workflow(1, 2, 3) == 9
    assert dynamic_workflow(1, 2, 3) == [6, 6, 6]


def test_decorators_v3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def add(a, b):
        return a + b

    def mult(a, b):
        return a * b

    def make_more(val):
        return [val] * 3

    def add_distributed(vals, c):
        return [job(add)(val, c) for val in vals]

    def workflow(a, b, c):
        return job(mult)(job(add)(a, b), c)

    def dynamic_workflow(a, b, c):
        result1 = job(add)(a, b)
        result2 = job(make_more)(result1)
        return subflow(add_distributed)(result2, c)

    assert job(add)(1, 2) == 3
    assert job(mult)(1, 2) == 2
    assert flow(workflow)(1, 2, 3) == 9
    assert flow(dynamic_workflow)(1, 2, 3) == [6, 6, 6]
