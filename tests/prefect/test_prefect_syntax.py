import pytest

from quacc import SETTINGS, flow, job, subflow

prefect = pytest.importorskip("prefect")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "prefect",
    reason="This test requires the Prefect workflow engine",
)
from prefect.testing.utilities import prefect_test_harness


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


def test_patch():
    from prefect import flow

    from quacc import job

    @job
    def task1(val: float) -> dict:
        return {"input": val, "result": val * 100}

    @job
    def task2(val: float) -> dict:
        return {"input": val, "result": val * 200}

    @flow
    def workflow(val):
        future1 = task1.submit(val)
        return task2(future1["result"])

    assert workflow(1) == {"input": 100, "result": 20000}


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

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2).result() == 3
    assert workflow(1, 2, 3).result() == 9
    results = dynamic_workflow(1, 2, 3)
    assert [result.result() for result in results] == [6, 6, 6]
