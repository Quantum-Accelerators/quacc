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
    if SETTINGS.WORKFLOW_ENGINE == "dask":
        func = strip_decorator(func)
        updated_func = job(partial(func, **params))
    else:
        updated_func = partial(func, **params)

    return updated_func


def customize_funcs(
    names: list[str] | str,
    funcs: list[Callable] | Callable,
    parameters: dict[str, dict[str, Any]] | None = None,
    decorators: dict[str, Callable | None] | None = None,
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

    Returns
    -------
    tuple[Callable] | Callable
        The customized functions, returned in the same order as provided in `funcs`.
    """
    parameters = parameters or {}
    decorators = decorators or {}
    updated_funcs = []

    if not isinstance(names, (list, tuple)):
        names = [names]
    if not isinstance(funcs, (list, tuple)):
        funcs = [funcs]

    funcs_dict = dict(zip(names, funcs))

    if "all" in names:
        raise ValueError("Invalid function name: 'all' is a reserved name.")
    if bad_decorator_keys := [k for k in decorators if k not in names and k != "all"]:
        raise ValueError(
            f"Invalid decorator keys: {bad_decorator_keys}. " f"Valid keys are: {names}"
        )
    if bad_parameter_keys := [k for k in parameters if k not in names and k != "all"]:
        raise ValueError(
            f"Invalid parameter keys: {bad_parameter_keys}. " f"Valid keys are: {names}"
        )

    for func_name, func in funcs_dict.items():
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

    if len(updated_funcs) == 1:
        return updated_funcs[0]

    return tuple(updated_funcs)
