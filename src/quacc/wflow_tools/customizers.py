"""Functions to customize workflow steps."""
from __future__ import annotations

from copy import deepcopy
from functools import partial
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, Callable


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
    from quacc import SETTINGS

    if SETTINGS.WORKFLOW_ENGINE == "covalent":
        if hasattr(func, "electron_object"):
            func = func.electron_object.function

        if hasattr(func, "workflow_function"):
            func = func.workflow_function.get_deserialized()

    elif SETTINGS.WORKFLOW_ENGINE == "dask":
        if hasattr(func, "__wrapped__"):
            func = func.__wrapped__

    elif SETTINGS.WORKFLOW_ENGINE == "jobflow":
        if hasattr(func, "original"):
            func = func.original

    elif SETTINGS.WORKFLOW_ENGINE in ("parsl", "redun"):
        if hasattr(func, "func"):
            func = func.func

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


def update_parameters(func: Callable, params: dict[str, Any]) -> Callable:
    """
    Update the parameters of a function.

    Parameters
    ----------
    func
        The function to update.
    params
        The parameters and associated values to update.

    Returns
    -------
    Callable
        The updated function.
    """
    from quacc import SETTINGS, job

    if SETTINGS.WORKFLOW_ENGINE != "dask":
        return partial(func, **params)

    func = strip_decorator(func)
    return job(partial(func, **params))


def customize_funcs(
    names: list[str] | str,
    funcs: list[Callable] | Callable,
    parameters: dict[str, dict[str, Any]] | None = None,
    decorators: dict[str, Callable | None] | None = None,
    prep_for_subflow: bool = True,
) -> tuple[Callable] | Callable:
    """
    Customize a set of functions with decorators and common parameters.

    Parameters
    ----------
    names
        The names of the functions to customize, in the order they should be returned.
    funcs
        The functions to customize, in the order they are described in `names`.
    parameters
        Custom parameters to apply to each function. The keys of this dictionary correspond
        to the strings in `names`. If the key `"all"` is present, it will be applied to all
        functions. If the value is `None`, no custom parameters will be applied to that function.
    decorators
        Custom decorators to apply to each function. The keys of this dictionary correspond
        to the strings in `names`. If the key `"all"` is present, it will be applied to all
        functions. If a value is `None`, no decorator will be applied that function.
    prep_for_subflow
        Whether to prepare the functions for use in a subflow. This is used to handle
        some peculiarities of the various workflow engines.

    Returns
    -------
    tuple[Callable] | Callable
        The customized functions, returned in the same order as provided in `funcs`.
    """
    from quacc import SETTINGS

    parameters = parameters or {}
    decorators = decorators or {}
    updated_funcs = []

    if not isinstance(names, (list, tuple)):
        names = [names]
    if not isinstance(funcs, (list, tuple)):
        funcs = [funcs]

    if "all" in names:
        raise ValueError("Invalid function name: 'all' is a reserved name.")
    if bad_decorator_keys := [k for k in decorators if k not in names and k != "all"]:
        raise ValueError(
            f"Invalid decorator keys: {bad_decorator_keys}. Valid keys are: {names}"
        )
    if bad_parameter_keys := [k for k in parameters if k not in names and k != "all"]:
        raise ValueError(
            f"Invalid parameter keys: {bad_parameter_keys}. Valid keys are: {names}"
        )

    for i, func in enumerate(funcs):
        func_ = deepcopy(func)
        if params := parameters.get("all"):
            func_ = update_parameters(func_, params)
        if params := parameters.get(names[i]):
            func_ = update_parameters(func_, params)
        if "all" in decorators:
            func_ = redecorate(func_, decorators["all"])
        if names[i] in decorators:
            func_ = redecorate(func_, decorators[names[i]])
        updated_funcs.append(func_)

    if prep_for_subflow and SETTINGS.WORKFLOW_ENGINE == "dask":
        from dask.core import literal
        from dask.delayed import Delayed

        updated_funcs = [
            literal(f) if isinstance(f, Delayed) else f for f in updated_funcs
        ]

    return updated_funcs[0] if len(updated_funcs) == 1 else tuple(updated_funcs)
