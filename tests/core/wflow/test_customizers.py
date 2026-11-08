from __future__ import annotations

from types import SimpleNamespace

import pytest

from quacc import Remove, job
from quacc.wflow_tools.customizers import (
    customize_jobs,
    resolve_job_params,
    update_parameters,
)


def test_customize_jobs():
    def add(a, b=1):
        return a + b

    def mult(a, b=2):
        return a * b

    jobs = customize_jobs(
        {"add": add, "mult": mult}, param_swaps={"all": {"b": 3}, "add": {"b": 4}}
    )
    assert list(jobs) == ["add", "mult"]
    assert jobs["add"](1) == 5
    assert jobs["mult"](2) == 6

    with pytest.raises(ValueError, match="Invalid parameter keys"):
        customize_jobs({"add": add}, param_swaps={"bad": {"b": 2}})

    with pytest.raises(ValueError, match="reserved name"):
        customize_jobs({"all": add})


def test_basic_customizers():
    def add(a, b=1, c=2):
        return a + b + c

    def mult(a, b=1, c=2, d=2):
        return a * b * c * d

    jobs = customize_jobs({"add": add, "mult": mult}, param_swaps={"add": {"b": 2}})
    add_, mult_ = jobs.values()
    assert add_(1) == 5
    assert mult_(1) == 4
    assert add(1) == 4
    assert mult(1) == 4

    jobs = customize_jobs(
        {"add": add, "mult": mult}, param_swaps={"add": {"b": 2}, "mult": {"b": 2}}
    )
    add_, mult_ = jobs.values()
    assert add_(1) == 5
    assert mult_(1) == 8
    assert add(1) == 4
    assert mult(1) == 4

    jobs = customize_jobs(
        {"add": add, "mult": mult},
        param_swaps={"add": {"b": 2}, "mult": {"b": 2}},
        decorators={"add": job(), "mult": job()},
    )
    add_, mult_ = jobs.values()
    assert add_(1) == 5
    assert mult_(1) == 8
    assert add(1) == 4
    assert mult(1) == 4

    jobs = customize_jobs({"add": add, "mult": mult}, param_swaps={"all": {"b": 2}})
    add_, mult_ = jobs.values()
    assert add_(1) == 5
    assert mult_(1) == 8
    assert add(1) == 4
    assert mult(1) == 4

    with pytest.raises(
        ValueError,
        match=r"Invalid parameter keys: \['bad'\]. Valid keys are: \['add', 'mult']",
    ):
        customize_jobs({"add": add, "mult": mult}, param_swaps={"bad": {"b": 2}})

    with pytest.raises(
        ValueError, match="Invalid function name: 'all' is a reserved name"
    ):
        customize_jobs({"all": add}, param_swaps={"all": {"b": 2}})

    with pytest.raises(ValueError, match="Invalid decorator keys"):
        customize_jobs({"add": add}, decorators={"bad": job})

    assert customize_jobs({"add": add})["add"](1) == 4


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


def test_resolve_job_params_precedence():
    params = resolve_job_params(
        ["relax_job", "static_job"],
        param_defaults={
            "all": {"a": "default_all"},
            "relax_job": {"a": "default_relax", "b": "default_relax"},
        },
        param_swaps={"all": {"a": "user_all"}, "static_job": {"a": "user_static"}},
    )
    # User "all" beats a recipe's job-specific default
    assert params["relax_job"] == {"a": "user_all", "b": "default_relax"}
    assert params["static_job"] == {"a": "user_static"}

    assert resolve_job_params(["relax_job"]) == {"relax_job": {}}

    with pytest.raises(ValueError, match="Invalid parameter keys"):
        resolve_job_params(["relax_job"], param_defaults={"bad": {}})


def test_customize_jobs_user_all_beats_default():
    def relax(atoms, relax_cell=False):
        return relax_cell

    jobs = customize_jobs(
        {"relax_job": relax},
        param_defaults={"relax_job": {"relax_cell": True}},
        param_swaps={"all": {"relax_cell": False}},
    )
    assert jobs["relax_job"](None) is False


def test_customize_jobs_nested_merge():
    def relax(atoms, opt_params=None):
        return opt_params

    jobs = customize_jobs(
        {"relax_job": relax},
        param_defaults={"relax_job": {"opt_params": {"fmax": 0.1, "max_steps": 100}}},
        param_swaps={
            "all": {"opt_params": {"fmax": 1e-3}},
            "relax_job": {"opt_params": {"max_steps": 500}},
        },
    )
    assert jobs["relax_job"](None) == {"fmax": 1e-3, "max_steps": 500}


def test_customize_jobs_remove_default():
    def relax(atoms, opt_params=None):
        return opt_params

    jobs = customize_jobs(
        {"relax_job": relax},
        param_defaults={"relax_job": {"opt_params": {"fmax": 0.1, "max_steps": 100}}},
        param_swaps={"relax_job": {"opt_params": {"max_steps": Remove}}},
    )
    assert jobs["relax_job"](None) == {"fmax": 0.1}


def test_customize_jobs_signature_check():
    def relax(atoms, fmax=0.1):
        return fmax

    def static(atoms):
        return atoms

    def flexible(atoms, **kwargs):
        return kwargs

    with pytest.raises(ValueError, match=r"Invalid parameters for 'relax_job': \['fmx'\]"):
        customize_jobs({"relax_job": relax}, param_swaps={"relax_job": {"fmx": 0.01}})

    with pytest.raises(ValueError, match=r"Invalid parameters for 'relax_job'"):
        customize_jobs({"relax_job": relax}, param_defaults={"relax_job": {"fmx": 0.01}})

    # "all" may contain parameters that only some jobs accept
    jobs = customize_jobs(
        {"relax_job": relax, "static_job": static}, param_swaps={"all": {"fmax": 0.01}}
    )
    assert jobs["relax_job"](None) == 0.01

    # **kwargs functions are not checked
    jobs = customize_jobs({"job": flexible}, param_swaps={"job": {"anything": 1}})
    assert jobs["job"](None) == {"anything": 1}


def test_customize_jobs_decorator_precedence():
    def add(a, b=1):
        return a + b

    calls = []

    def tag(label):
        def decorator(func):
            def wrapper(*args, **kwargs):
                calls.append(label)
                return func(*args, **kwargs)

            return wrapper

        return decorator

    jobs = customize_jobs(
        {"add": add, "add2": add},
        param_swaps={"all": {"b": 2}},
        decorators={"all": tag("all"), "add": tag("specific")},
    )
    assert jobs["add"](1) == 3
    assert jobs["add2"](1) == 3
    assert calls == ["specific", "all"]
