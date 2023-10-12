import pytest

from quacc import SETTINGS

redun = pytest.importorskip("redun")

DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    SETTINGS.WORKFLOW_ENGINE = "redun"


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


@pytest.fixture()
def scheduler():
    return redun.Scheduler()


def test_redun_decorators(tmpdir, scheduler):
    tmpdir.chdir()
    from quacc import flow, job, subflow

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

    assert scheduler.run(add(1, 2)) == 3
    assert scheduler.run(mult(1, 2)) == 2
    assert scheduler.run(workflow(1, 2, 3)) == 9
    assert scheduler.run(dynamic_workflow(1, 2, 3)) == [6, 6, 6]
