from __future__ import annotations

from types import SimpleNamespace

import pytest

from quacc import job
from quacc.wflow_tools.customizers import customize_funcs, update_parameters


def test_basic_customizers():
    def add(a, b=1, c=2):
        return a + b + c

    def mult(a, b=1, c=2, d=2):
        return a * b * c * d

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], param_swaps={"add": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 4
    assert add(1) == 4
    assert mult(1) == 4

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], param_swaps={"add": {"b": 2}, "mult": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 8
    assert add(1) == 4
    assert mult(1) == 4

    add_, mult_ = customize_funcs(
        ["add", "mult"],
        [add, mult],
        param_swaps={"add": {"b": 2}, "mult": {"b": 2}},
        decorators={"add": job(), "mult": job()},
    )
    assert add_(1) == 5
    assert mult_(1) == 8
    assert add(1) == 4
    assert mult(1) == 4

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], param_swaps={"all": {"b": 2}}
    )
    assert add_(1) == 5
    assert mult_(1) == 8
    assert add(1) == 4
    assert mult(1) == 4

    with pytest.raises(
        ValueError,
        match=r"Invalid parameter keys: \['bad'\]. Valid keys are: \['add', 'mult']",
    ):
        customize_funcs(["add", "mult"], [add, mult], param_swaps={"bad": {"b": 2}})

    with pytest.raises(
        ValueError, match="Invalid function name: 'all' is a reserved name"
    ):
        customize_funcs("all", [add], param_swaps={"all": {"b": 2}})

    with pytest.raises(ValueError, match="Invalid decorator keys"):
        customize_funcs("add", add, decorators={"bad": job})

    assert customize_funcs("add", add)(1) == 4


def test_update_parameters_edge_cases(monkeypatch):
    class CallableWithoutName:
        def __call__(self, value=1):
            return value

    updated = update_parameters(CallableWithoutName(), {}, decorator=None)
    assert updated.__name__ == ""
    assert updated() == 1

    monkeypatch.setattr(
        "quacc.get_settings", lambda: SimpleNamespace(WORKFLOW_ENGINE="dask")
    )
    with pytest.raises(ValueError, match="Invalid decorator name: invalid"):
        update_parameters(lambda: None, {}, decorator="invalid")
