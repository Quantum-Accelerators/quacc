"""Functions to customize workflow steps."""
from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

from quacc import SETTINGS, job

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

    if hasattr(func, "__wrapped__"):
        func = func.__wrapped__

    if SETTINGS.WORKFLOW_ENGINE == "covalent":
        from covalent._workflow.lattice import Lattice

        if isinstance(func, Lattice):
            func = func.workflow_function.get_deserialized()
    return func


def redecorate(func: Callable, decorator: Callable | None) -> Callable:
    """
    Redecorate pre-decorated functions with custom decorators.

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
    Update the parameters of a function. If the function does not have a given parameter,
    it is ignored.

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
    if SETTINGS.WORKFLOW_ENGINE == "dask":
        func = strip_decorator(func)
        updated_func = job(partial(func, **params))
    else:
        updated_func = partial(func, **params)

    return updated_func


def customize_funcs(
    funcs: dict[str, Callable],
    decorators: dict[str, Callable | None] | None = None,
    parameters: dict[str, Any] | None = None,
) -> tuple[Callable]:
    """
    Customize a set of functions with decorators and common parameters.

    Parameters
    ----------
    funcs
        The functions to customize, as a dictionary where the keys are unique
        identifiers for each function and the values are the functions themselves.
    decorators
        Custom decorators to apply to each function. The keys of this dictionary correspond
        to the keys of `funcs`. If the key `"all"` is present, it will be applied to all
        functions. If a value is `None`, no decorator will be applied that function.
    parameters
        Custom parameters to apply to each function. The keys of this dictionary correspond
        to the keys of `funcs`. If the key `"all"` is present, it will be applied to all
        functions. If the value is `None`, no custom parameters will be applied to that function.

    Returns
    -------
    tuple[Callable]
        The customized functions, returned in the same order as provided in `funcs`.
    """
    decorators = decorators or {}
    parameters = parameters or {}
    updated_funcs = []

    if bad_decorator_keys := [k for k in decorators if k not in funcs and k != "all"]:
        raise ValueError(
            f"Invalid decorator keys: {bad_decorator_keys}. "
            f"Valid keys are: {list(funcs.keys())}"
        )
    if bad_parameter_keys := [k for k in parameters if k not in funcs and k != "all"]:
        raise ValueError(
            f"Invalid parameter keys: {bad_parameter_keys}. "
            f"Valid keys are: {list(funcs.keys())}"
        )

    for func_name, func in funcs.items():
        func_ = func
        if params := parameters.get("all"):
            func_ = update_parameters(func_, params)
        if params := parameters.get(func_name):
            func_ = update_parameters(func_, params)
        if "all" in decorators:
            func_ = redecorate(func_, decorators["all"])
        if func_name in decorators:
            func_ = redecorate(func_, decorators[func_name])
        updated_funcs.append(func_)
    return tuple(updated_funcs)
