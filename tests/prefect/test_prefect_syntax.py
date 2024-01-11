import pytest

prefect = pytest.importorskip("prefect")

from quacc import flow, job, strip_decorator, subflow


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

    assert workflow(1).result() == {"input": 100, "result": 20000}


def test_patch_local():
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    SETTINGS.PREFECT_AUTO_SUBMIT = False

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
    SETTINGS.PREFECT_AUTO_SUBMIT = DEFAULT_SETTINGS.PREFECT_AUTO_SUBMIT


def test_prefect_decorators(tmp_path, monkeypatch):
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

    @flow
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2).is_completed()
    assert workflow(1, 2, 3).is_completed()
    assert [r.result() for r in dynamic_workflow(1, 2, 3)] == [6, 6, 6]
    assert [r.result() for r in dynamic_workflow2(1, 2, 3)] == [6, 6, 6]
    assert dynamic_workflow3(1, 2, 3).result() == 12


def test_prefect_decorators_local(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    SETTINGS.PREFECT_AUTO_SUBMIT = False

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

    SETTINGS.PREFECT_AUTO_SUBMIT = DEFAULT_SETTINGS.PREFECT_AUTO_SUBMIT


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
