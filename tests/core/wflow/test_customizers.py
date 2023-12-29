import pytest

from quacc.wflow_tools.customizers import customize_funcs


def test_basic_customizers():
    def add(a, b=1, c=2):
        return a + b + c

    def mult(a, b=1, c=2, d=2):
        return a * b * c * d

    add_, mult_ = customize_funcs(
        {"add": add, "mult": mult}, parameters={"add": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 4

    add_, mult_ = customize_funcs(
        {"add": add, "mult": mult}, parameters={"add": {"b": 2}, "mult": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 8

    add_, mult_ = customize_funcs(
        {"add": add, "mult": mult}, parameters={"all": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 8

    with pytest.raises(ValueError):
        customize_funcs({"add": add, "mult": mult}, parameters={"bad": {"b": 2}})
