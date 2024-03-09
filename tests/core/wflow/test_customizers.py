import pytest

from quacc import flow, job
from quacc.wflow_tools.customizers import customize_funcs


def test_basic_customizers():
    def add(a, b=1, c=2):
        return a + b + c

    def mult(a, b=1, c=2, d=2):
        return a * b * c * d

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], parameters={"add": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 4

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], parameters={"add": {"b": 2}, "mult": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 8

    add_, mult_ = customize_funcs(
        ["add", "mult"],
        [add, mult],
        parameters={"add": {"b": 2}, "mult": {"b": 2}},
        decorators={"add": job(), "mult": job()},
    )
    assert add_(1) == 5
    assert mult_(1) == 8

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], parameters={"all": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 8

    with pytest.raises(ValueError):
        customize_funcs(["add", "mult"], [add, mult], parameters={"bad": {"b": 2}})

    with pytest.raises(ValueError):
        customize_funcs("all", [add], parameters={"all": {"b": 2}})
