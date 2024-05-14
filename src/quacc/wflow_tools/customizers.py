"""Functions to customize workflow steps."""

from __future__ import annotations

from copy import deepcopy
from functools import partial
from typing import TYPE_CHECKING

from quacc.utils.dicts import recursive_dict_merge

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

    decorator = getattr(func, "quacc_decorator", None)
    if decorator is None:
        return func

    if SETTINGS.WORKFLOW_ENGINE == "covalent":
        if decorator in ("job", "subflow"):
            func = func.electron_object.function

        if decorator in ("flow", "subflow"):
            func = func.workflow_function.get_deserialized()

    elif SETTINGS.WORKFLOW_ENGINE == "dask":
        if decorator == "job":
            func = func.func
        func = func.__wrapped__
        if decorator == "subflow":
            func = func.__wrapped__

    elif SETTINGS.WORKFLOW_ENGINE == "jobflow":
        func = func.original

    elif SETTINGS.WORKFLOW_ENGINE == "parsl":
        func = func.func

    elif SETTINGS.WORKFLOW_ENGINE == "prefect":
        if SETTINGS.PREFECT_AUTO_SUBMIT:
            func = func.__wrapped__
        func = func.fn

    elif SETTINGS.WORKFLOW_ENGINE == "redun":
        func = func.func

    return func


def redecorate(func: Callable, decorator: Callable) -> Callable:
    """
    Redecorate a pre-decorated function with a custom decorator.

    Parameters
    ----------
    func
        The pre-decorated function.
    decorator
        The new decorator to apply.

    Returns
    -------
    Callable
        The newly decorated function.
    """
    func = strip_decorator(func)
    return decorator(func)


def update_parameters(func: Callable, params: dict[str, Any]) -> Callable:
    """
    Update the parameters of a (potentially decorated) function.

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
    from quacc import SETTINGS, flow, job, subflow

    if (
        decorators_type := hasattr(func, "quacc_decorator")
        and SETTINGS.WORKFLOW_ENGINE == "dask"
    ):
        if decorators_type == "job":
            decorator = job
        elif decorators_type == "flow":
            decorator = flow
        elif decorators_type == "subflow":
            decorator = subflow

        func = strip_decorator(func)
        return decorator(partial(func, **params))

    return partial(func, **params)


def customize_funcs(
    names: list[str] | str,
    funcs: list[Callable] | Callable,
    param_defaults: dict[str, dict[str, Any]] | None = None,
    param_swaps: dict[str, dict[str, Any]] | None = None,
    decorators: dict[str, Callable | None] | None = None,
) -> tuple[Callable, ...] | Callable:
    """
    Customize a set of functions with decorators and common parameters.

    Parameters
    ----------
    names
        The names of the functions to customize, in the order they should be returned.
    funcs
        The functions to customize, in the order they are described in `names`.
    param_defaults
        Default parameters to apply to each function. The keys of this dictionary correspond
        to the strings in `names`. If the key `"all"` is present, it will be applied to all
        functions. If the value is `None`, no custom parameters will be applied to that function.
    param_swaps
        User-overrides of parameters to apply to each function. The keys of this dictionary correspond
        to the strings in `names`. If the key `"all"` is present, it will be applied to all
        functions. If the value is `None`, no custom parameters will be applied to that function.
    decorators
        Custom decorators to apply to each function. The keys of this dictionary correspond
        to the strings in `names`. If the key `"all"` is present, it will be applied to all
        functions. If a value is `None`, no decorator will be applied that function.

    Returns
    -------
    tuple[Callable, ...] | Callable
        The customized functions, returned in the same order as provided in `funcs`.
    """
    parameters = recursive_dict_merge(param_defaults, param_swaps)
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
        if decorator := decorators.get("all"):
            func_ = redecorate(func_, decorator)
        if decorator := decorators.get(names[i]):
            func_ = redecorate(func_, decorator)
        if params := parameters.get("all"):
            func_ = update_parameters(func_, params)
        if params := parameters.get(names[i]):
            func_ = update_parameters(func_, params)
        updated_funcs.append(func_)

    return updated_funcs[0] if len(updated_funcs) == 1 else tuple(updated_funcs)
