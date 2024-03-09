"""Functions to customize workflow steps."""

from __future__ import annotations

from copy import deepcopy
from functools import partial
from typing import TYPE_CHECKING, Literal

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
        from covalent._workflow.lattice import Lattice

        if hasattr(func, "electron_object"):
            func = func.electron_object.function

        if isinstance(func, Lattice):
            func = func.workflow_function.get_deserialized()

    elif SETTINGS.WORKFLOW_ENGINE == "dask":
        from dask.delayed import Delayed

        from quacc.wflow_tools.decorators import Delayed_

        if isinstance(func, Delayed_):
            func = func.func
        if isinstance(func, Delayed):
            func = func.__wrapped__
            if hasattr(func, "__wrapped__"):
                # Needed for custom `@subflow` decorator
                func = func.__wrapped__

    elif SETTINGS.WORKFLOW_ENGINE == "jobflow":
        if hasattr(func, "original"):
            func = func.original

    elif SETTINGS.WORKFLOW_ENGINE == "parsl":
        from parsl.app.python import PythonApp

        if isinstance(func, PythonApp):
            func = func.func

    elif SETTINGS.WORKFLOW_ENGINE == "prefect":
        from prefect import Flow as PrefectFlow
        from prefect import Task

        if isinstance(func, (Task, PrefectFlow)):
            func = func.fn
        elif hasattr(func, "__wrapped__"):
            func = func.__wrapped__

    elif SETTINGS.WORKFLOW_ENGINE == "redun":
        from redun import Task

        if isinstance(func, Task):
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
    from quacc import SETTINGS, flow, job, subflow

    if decorator and SETTINGS.WORKFLOW_ENGINE == "dask":
        if decorator == "job":
            decorator = job
        elif decorator == "flow":
            decorator = flow
        elif decorator == "subflow":
            decorator = subflow

        func = strip_decorator(func)
        return decorator(partial(func, **params))

    return partial(func, **params)


def customize_funcs(
    names: list[str] | str,
    funcs: list[Callable] | Callable,
    parameters: dict[str, dict[str, Any]] | None = None,
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
    tuple[Callable, ...] | Callable
        The customized functions, returned in the same order as provided in `funcs`.
    """

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
        if "all" in decorators:
            func_ = redecorate(func_, decorators["all"])
        if names[i] in decorators:
            func_ = redecorate(func_, decorators[names[i]])
        if params := parameters.get("all"):
            func_ = update_parameters(func_, params)
        if params := parameters.get(names[i]):
            func_ = update_parameters(func_, params)
        updated_funcs.append(func_)

    return updated_funcs[0] if len(updated_funcs) == 1 else tuple(updated_funcs)
