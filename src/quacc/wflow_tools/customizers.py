"""Functions to customize workflow steps."""
from __future__ import annotations

import inspect
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
    if decorator is None:
        return func
    return decorator(func)


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
    stripped_func = strip_decorator(func)
    func_params = inspect.signature(stripped_func).parameters
    valid_params = {k: v for k, v in params.items() if k in func_params}
    return partial(func, **valid_params)


def customize_funcs(
    funcs: dict[str, Callable],
    decorators: dict[str, Callable | None] | None,
    parameters: dict[str, Any] | None,
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
        If a function does not have a given parameter, it is ignored.

    Returns
    -------
    tuple[Callable]
        The customized functions, returned in the same order as provided in `funcs`.
    """
    decorators = decorators or {}
    parameters = parameters or {}
    updated_funcs = []
    for func_name, func in funcs.items():
        if params := parameters.get("all"):
            func = update_parameters(func, params)
        if params := parameters.get(func_name):
            func = update_parameters(func, params)
        if "all" in decorators:
            func = redecorate(func, decorators["all"])
        if func_name in decorators:
            func = redecorate(func, decorators[func_name])
        updated_funcs.append(func)
    return tuple(updated_funcs)
