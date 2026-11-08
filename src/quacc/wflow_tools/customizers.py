"""Functions to customize workflow steps."""

from __future__ import annotations

import inspect
from copy import deepcopy
from functools import partial
from typing import TYPE_CHECKING, Literal

from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
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
        from prefect import Task as PrefectTask

        if isinstance(func, PrefectTask | PrefectFlow):
            func = func.fn
        elif hasattr(func, "__wrapped__"):
            func = func.__wrapped__

    elif settings.WORKFLOW_ENGINE == "redun":
        from redun import Task as RedunTask

        if isinstance(func, RedunTask):
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
        vars(partial_fn)["__name__"] = func.name
    elif hasattr(func, "__name__"):
        vars(partial_fn)["__name__"] = func.__name__
    else:
        vars(partial_fn)["__name__"] = ""

    return partial_fn


def _check_keys(kind: str, given: Iterable[str], names: list[str]) -> None:
    """Raise if `given` contains keys that are neither a job name nor "all"."""
    if bad_keys := [k for k in given if k not in names and k != "all"]:
        raise ValueError(f"Invalid {kind} keys: {bad_keys}. Valid keys are: {names}")


def resolve_job_params(
    names: Iterable[str],
    param_defaults: dict[str, dict[str, Any] | None] | None = None,
    param_swaps: dict[str, dict[str, Any] | None] | None = None,
) -> dict[str, dict[str, Any]]:
    """
    Resolve the final keyword arguments for each named job in a workflow.

    The precedence is (lowest to highest): `param_defaults["all"]`,
    `param_defaults[name]`, `param_swaps["all"]`, `param_swaps[name]`. Nested
    dictionaries are merged recursively at every step, and any value set to
    [quacc.Remove][] is dropped from the result.

    Parameters
    ----------
    names
        The job names used by the workflow.
    param_defaults
        Default (i.e. recipe-provided) parameters, keyed by job name or `"all"`.
    param_swaps
        User overrides (i.e. `job_params`), keyed by job name or `"all"`.

    Returns
    -------
    dict[str, dict[str, Any]]
        The fully merged parameters for every name in `names`.
    """
    names = list(names)
    if "all" in names:
        raise ValueError("Invalid function name: 'all' is a reserved name.")
    param_defaults = param_defaults or {}
    param_swaps = param_swaps or {}
    _check_keys("parameter", param_defaults, names)
    _check_keys("parameter", param_swaps, names)

    return {
        name: recursive_dict_merge(
            param_defaults.get("all"),
            param_defaults.get(name),
            param_swaps.get("all"),
            param_swaps.get(name),
        )
        for name in names
    }


def _check_signature(name: str, func: Callable, keys: Iterable[str]) -> None:
    """
    Raise if `func` cannot accept the keyword arguments `keys`.

    This catches typos when the workflow is built rather than when the job runs.
    Functions that accept `**kwargs`, or that cannot be introspected, are skipped.
    """
    try:
        signature = inspect.signature(strip_decorator(func))
    except (TypeError, ValueError):
        return
    if any(p.kind is p.VAR_KEYWORD for p in signature.parameters.values()):
        return
    if unknown := sorted(set(keys) - set(signature.parameters)):
        raise ValueError(
            f"Invalid parameters for '{name}': {unknown}. "
            f"Valid parameters are: {list(signature.parameters)}"
        )


def customize_jobs(
    jobs: dict[str, Callable],
    param_defaults: dict[str, dict[str, Any] | None] | None = None,
    param_swaps: dict[str, dict[str, Any] | None] | None = None,
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
        Default parameters to apply to each job, typically set by the recipe.
        The keys correspond to the keys in `jobs`, or `"all"` to apply to every
        job.
    param_swaps
        User overrides of parameters to apply to each job, typically the
        recipe's `job_params`. The keys correspond to the keys in `jobs`, or
        `"all"` to apply to every job. See [quacc.wflow_tools.customizers.resolve_job_params][]
        for how these are merged with `param_defaults`.
    decorators
        Custom decorators to apply to each job. The keys correspond to the keys
        in `jobs`, or `"all"` to apply to every job. A job-specific decorator
        takes precedence over `"all"`. A value of `None` strips the existing
        decorator from that job.

    Returns
    -------
    dict[str, Callable]
        Mapping of customized jobs in the same order as `jobs`.
    """
    names = list(jobs)
    decorators = decorators or {}
    parameters = resolve_job_params(names, param_defaults, param_swaps)
    _check_keys("decorator", decorators, names)

    updated_jobs = {}
    for name, func in jobs.items():
        func_ = deepcopy(func)

        # Only job-specific keys are checked: "all" may hold parameters that
        # some jobs in the workflow don't accept but never run with.
        job_specific_keys = {
            *((param_defaults or {}).get(name) or {}),
            *((param_swaps or {}).get(name) or {}),
        }
        _check_signature(name, func_, job_specific_keys)

        stripped = False
        if name in decorators or "all" in decorators:
            new_decorator = (
                decorators[name] if name in decorators else decorators["all"]
            )
            func_ = redecorate(func_, new_decorator)
            stripped = new_decorator is None

        if params := parameters[name]:
            func_ = update_parameters(
                func_, params, decorator=None if stripped else "job"
            )
        updated_jobs[name] = func_

    return updated_jobs
