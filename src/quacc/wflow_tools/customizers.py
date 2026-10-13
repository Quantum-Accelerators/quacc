"""Functions to customize workflow steps."""

from __future__ import annotations

from copy import deepcopy
from functools import partial
from typing import TYPE_CHECKING, Literal

from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any


def strip_decorator(func: Callable) -> Callable:
    """
    Strip the decorators from a function.

    Parameters
    ----------
    func
        The function to strip decorators from.

    Returns
    -------
    Callable
        The function with all decorators removed.
    """
    from quacc import get_settings

    settings = get_settings()

    if settings.WORKFLOW_ENGINE == "dask":
        from dask.delayed import Delayed

        from quacc.wflow_tools.decorators import Delayed_

        if isinstance(func, Delayed_):
            func = func.func
        if isinstance(func, Delayed):
            func = func.__wrapped__
            if hasattr(func, "__wrapped__"):
                # Needed for custom `@subflow` decorator
                func = func.__wrapped__

    elif settings.WORKFLOW_ENGINE == "jobflow":
        if hasattr(func, "original"):
            func = func.original

    elif settings.WORKFLOW_ENGINE == "parsl":
        from parsl.app.python import PythonApp

        if isinstance(func, PythonApp):
            func = func.func

    elif settings.WORKFLOW_ENGINE == "prefect":
        from prefect import Flow as PrefectFlow
        from prefect import Task

        if isinstance(func, Task | PrefectFlow):
            func = func.fn
        elif hasattr(func, "__wrapped__"):
            func = func.__wrapped__

    elif settings.WORKFLOW_ENGINE == "redun":
        from redun import Task

        if isinstance(func, Task):
            func = func.func

    elif settings.WORKFLOW_ENGINE == "ray":
        if hasattr(func, "__wrapped__"):
            func = func.__wrapped__

    return func


def redecorate(func: Callable, decorator: Callable | None) -> Callable:
    """
    Redecorate a pre-decorated function with a custom decorator.

    Parameters
    ----------
    func
        The pre-decorated function.
    decorator
        The new decorator to apply. If `None`, the function is stripped of its
        decorators.

    Returns
    -------
    Callable
        The newly decorated function.
    """
    func = strip_decorator(func)
    return func if decorator is None else decorator(func)


def update_parameters(
    func: Callable,
    params: dict[str, Any],
    decorator: Literal["job", "flow", "subflow"] | None = "job",
) -> Callable:
    """
    Update the parameters of a (potentially decorated) function.

    Parameters
    ----------
    func
        The function to update.
    params
        The parameters and associated values to update.
    decorator
        The decorator associated with `func`.

    Returns
    -------
    Callable
        The updated function.
    """
    from quacc import flow, get_settings, job, subflow

    settings = get_settings()

    if decorator and settings.WORKFLOW_ENGINE in ("dask", "ray"):
        if decorator == "job":
            decorator_func = job
        elif decorator == "flow":
            decorator_func = flow
        elif decorator == "subflow":
            decorator_func = subflow
        else:
            raise ValueError(
                f"Invalid decorator name: {decorator}. Valid names are: 'job', 'flow', 'subflow'"
            )
        func = strip_decorator(func)
        return decorator_func(partial(func, **params))

    partial_fn = partial(func, **params)
    # Assigning a __name__ allows monty's jsanitize function to work correctly
    # with this partial function.
    if hasattr(func, "name"):
        partial_fn.__name__ = func.name
    elif hasattr(func, "__name__"):
        partial_fn.__name__ = func.__name__
    else:
        partial_fn.__name__ = ""

    return partial_fn


def _get_parameter_update_decorator(
    func: Callable,
) -> Literal["job", "flow", "subflow"]:
    """Get the quacc decorator to preserve when updating parameters."""
    from quacc import get_settings

    if get_settings().WORKFLOW_ENGINE == "dask":
        from dask.delayed import Delayed

        from quacc.wflow_tools.decorators import Delayed_

        if isinstance(func, Delayed_):
            return "job"
        if isinstance(func, Delayed):
            return "subflow"
        if hasattr(func, "original"):
            return "flow"

    return "job"


def customize_jobs(
    jobs: dict[str, Callable],
    param_defaults: dict[str, dict[str, Any]] | None = None,
    param_swaps: dict[str, dict[str, Any]] | None = None,
    decorators: dict[str, Callable | None] | None = None,
) -> dict[str, Callable]:
    """
    Customize a mapping of workflow jobs with decorators and common parameters.

    Parameters
    ----------
    jobs
        Mapping of job names to functions. The names are used as keys in
        `param_defaults`, `param_swaps`, and `decorators`.
    param_defaults
        Default parameters to apply to each job. The keys of this dictionary
        correspond to the keys in `jobs`. If the key `"all"` is present, its
        parameters are applied to every job before job-specific parameters.
        If the value is `None`, no default parameters are applied to that job.
    param_swaps
        User overrides of parameters to apply to each job. The keys of this
        dictionary correspond to the keys in `jobs`. If the key `"all"` is
        present, its parameters are applied to every job before job-specific
        parameters. If the value is `None`, no parameters are applied to that job.
    decorators
        Custom decorators to apply to each job. The keys of this dictionary
        correspond to the keys in `jobs`. If the key `"all"` is present, its
        decorator is applied to every job before job-specific decorators. A value
        of `None` strips the existing decorator from that job.

    Returns
    -------
    dict[str, Callable]
        Mapping of customized jobs in the same order as `jobs`.
    """
    parameters = recursive_dict_merge(param_defaults, param_swaps)
    decorators = decorators or {}
    names = list(jobs)

    if "all" in jobs:
        raise ValueError("Invalid function name: 'all' is a reserved name.")
    if bad_decorator_keys := [k for k in decorators if k not in jobs and k != "all"]:
        raise ValueError(
            f"Invalid decorator keys: {bad_decorator_keys}. Valid keys are: {names}"
        )
    if bad_parameter_keys := [k for k in parameters if k not in jobs and k != "all"]:
        raise ValueError(
            f"Invalid parameter keys: {bad_parameter_keys}. Valid keys are: {names}"
        )

    updated_jobs = {}
    for name, func in jobs.items():
        func_ = deepcopy(func)
        if "all" in decorators:
            func_ = redecorate(func_, decorators["all"])
        if name in decorators:
            func_ = redecorate(func_, decorators[name])
        if params := parameters.get("all"):
            func_ = update_parameters(
                func_, params, decorator=_get_parameter_update_decorator(func_)
            )
        if params := parameters.get(name):
            func_ = update_parameters(
                func_, params, decorator=_get_parameter_update_decorator(func_)
            )
        updated_jobs[name] = func_

    return updated_jobs
