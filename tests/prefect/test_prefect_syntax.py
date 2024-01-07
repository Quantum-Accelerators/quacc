prefect = pytest.importorskip("prefect")
import pytest
from prefect.testing.utilities import prefect_test_harness

from quacc import flow, job, subflow


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


def test_patch():
    from quacc import SETTINGS

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

    assert workflow(1).is_completed()


def test_patch_local():
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    SETTINGS.PREFECT_TASK_RUNNER = False

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
    SETTINGS.PREFECT_TASK_RUNNER = DEFAULT_SETTINGS.PREFECT_TASK_RUNNER


def test_prefect_decorators(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from quacc import SETTINGS

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
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2).is_completed()
    assert workflow(1, 2, 3).is_completed()
    assert [r.is_completed() for r in dynamic_workflow(1, 2, 3)]
    assert [r.is_completed() for r in dynamic_workflow2(1, 2, 3)]


def test_prefect_decorators_local(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    SETTINGS.PREFECT_TASK_RUNNER = False

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
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2) == 3
    assert workflow(1, 2, 3) == 9
    results = dynamic_workflow(1, 2, 3)
    assert results == [6, 6, 6]
    results = dynamic_workflow2(1, 2, 3)
    assert results == [6, 6, 6]

    SETTINGS.PREFECT_TASK_RUNNER = DEFAULT_SETTINGS.PREFECT_TASK_RUNNER
