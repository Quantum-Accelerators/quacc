from __future__ import annotations

import pytest

redun = pytest.importorskip("redun")

from quacc import flow, job, subflow


@pytest.fixture()
def scheduler():
    return redun.Scheduler()


def test_redun_decorators(tmp_path, monkeypatch, scheduler):
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

    assert scheduler.run(add(1, 2)) == 3
    assert scheduler.run(mult(1, 2)) == 2
    assert scheduler.run(workflow(1, 2, 3)) == 9
    assert scheduler.run(dynamic_workflow(1, 2, 3)) == [6, 6, 6]
    assert scheduler.run(dynamic_workflow2(1, 2, 3)) == [6, 6, 6]
    assert scheduler.run(dynamic_workflow3(1, 2, 3)) == 12
