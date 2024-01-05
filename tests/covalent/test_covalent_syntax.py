import pytest

from quacc import SETTINGS, flow, job, strip_decorator, subflow

ct = pytest.importorskip("covalent")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "covalent",
    reason="This test requires the Covalent workflow engine",
)
from covalent._workflow.lattice import Lattice


def test_covalent_decorators(tmp_path, monkeypatch):
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

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert ct.get_result(ct.dispatch(workflow)(1, 2, 3), wait=True).result == 9
    assert ct.get_result(ct.dispatch(dynamic_workflow)(1, 2, 3), wait=True).result == [
        6,
        6,
        6,
    ]
    assert ct.get_result(ct.dispatch(dynamic_workflow2)(1, 2, 3), wait=True).result == [
        6,
        6,
        6,
    ]
    assert ct.get_result(
        ct.dispatch(flow(add_distributed))([1, 1, 1], 2), wait=True
    ).result == [3, 3, 3]


def test_covalent_decorators_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

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

    assert hasattr(add, "electron_object")
    stripped_add = strip_decorator(add)
    assert stripped_add(1, 2) == 3
    assert not hasattr(stripped_add, "electron_object")

    assert isinstance(add2, Lattice)
    stripped_add2 = strip_decorator(add2)
    assert stripped_add2(1, 2) == 3
    assert not hasattr(stripped_add2, "electron_object")
    assert not isinstance(stripped_add2, Lattice)

    assert hasattr(add3, "electron_object")
    stripped_add3 = strip_decorator(add3)
    assert stripped_add3(1, 2) == 3
    assert not hasattr(stripped_add3, "electron_object")
    assert not isinstance(stripped_add3, Lattice)
