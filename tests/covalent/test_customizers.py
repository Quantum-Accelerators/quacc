import pytest

ct = pytest.importorskip("covalent")


from quacc import SETTINGS, flow, job, subflow
from quacc.wflow_tools.customizers import customize_funcs


def test_customize_jobs():
    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], decorators={"add": job(executor="local")}
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor is None

    add_, mult_ = customize_funcs(
        ["add", "mult"],
        [add, mult],
        decorators={"add": job(executor="local"), "mult": job(executor="local")},
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor == "local"

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], decorators={"all": job(executor="local")}
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor == "local"

    add_, mult_ = customize_funcs(
        ["add", "mult"],
        [add, mult],
        decorators={"all": job(executor="local"), "mult": job(executor="dask")},
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor == "dask"


def test_customize_flows():
    @job
    def add(a, b):
        return a + b

    @flow
    def workflow():
        return add(1, 2)

    workflow_ = customize_funcs(
        "workflow", workflow, decorators={"workflow": flow(executor="local")}
    )
    assert workflow_.metadata["executor"] == "local"


def test_customize_subflows():
    @job
    def add(a, b):
        return a + b

    @subflow
    def subworkflow():
        return add(1, 2)

    subworkflow_ = customize_funcs(
        "workflow", subworkflow, decorators={"workflow": subflow(executor="local")}
    )
    assert subworkflow_.electron_object.executor == "local"


def test_bad_customizers():
    @job
    def add(a, b):
        return a + b

    with pytest.raises(ValueError):
        customize_funcs("add", add, decorators={"bad": job(executor="local")})
