import pytest

from quacc import SETTINGS

ct = pytest.importorskip("covalent")


DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    SETTINGS.WORKFLOW_ENGINE = "covalent"


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


def test_covalent_decorators(tmpdir):
    from quacc import flow, job, subflow

    tmpdir.chdir()

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
    assert ct.get_result(ct.dispatch(workflow)(1, 2, 3), wait=True).result == 9
    assert ct.get_result(ct.dispatch(dynamic_workflow)(1, 2, 3), wait=True).result == [
        6,
        6,
        6,
    ]
    assert ct.get_result(
        ct.dispatch(flow(add_distributed))([1, 1, 1], 2), wait=True
    ).result == [
        3,
        3,
        3,
    ]


def test_covalent_decorators_args(tmpdir):
    from quacc import flow, job, subflow

    tmpdir.chdir()

    @job(executor="local")
    def add(a, b):
        return a + b

    @job(executor="local")
    def mult(a, b):
        return a * b

    @job(executor="local")
    def make_more(val):
        return [val] * 3

    @subflow(executor="local")
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow(executor="local")
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow(executor="local")
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert ct.get_result(ct.dispatch(workflow)(1, 2, 3), wait=True).result == 9
    assert ct.get_result(ct.dispatch(dynamic_workflow)(1, 2, 3), wait=True).result == [
        6,
        6,
        6,
    ]
