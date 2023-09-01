from __future__ import annotations

import functools
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, Literal, TypeVar

    from ase import Atoms

    Job = TypeVar("Job")
    Flow = TypeVar("Flow")
    Subflow = TypeVar("Subflow")


def job(
    _func: callable | None = None,
    engine: Literal["covalent", "parsl", "jobflow", "local"] | None = None,
    **kwargs,
) -> Job:  # sourcery skip
    """
    Decorator for individual compute jobs. This is a @job decorator.

    job = ct.electron [Covalent] = python_app [Parsl] = job [Jobflow].

    The wrapped function gets a new kwarg, `decorator_kwargs`, that can
    be used to modify the workflow engine decorator keyword arguments.
    The wrapped function can also be accessed without the decorator
    by calling the `.__wrapped__` attribute.

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    engine
        The workflow engine. This defaults to `SETTINGS.WORKFLOW_ENGINE` if `None`.
    **kwargs
        Keyword arguments to pass to the workflow engine decorator.

    Returns
    -------
    Job
        The @job-decorated function.
    """

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        """
        The @job-decorated function.

        Parameters
        ----------
        *f_args
            Positional arguments to the function, if any.
        decorator_kwargs
            Keyword arguments to pass to the workflow engine decorator.
        **f_kwargs
            Keyword arguments to the function, if any.

        Returns
        -------
        Any
            The output of the @job-decorated function.
        """

        from quacc import SETTINGS

        if not decorator_kwargs:
            decorator_kwargs = kwargs

        wflow_engine = engine or SETTINGS.WORKFLOW_ENGINE
        if wflow_engine == "covalent":
            import covalent as ct

            decorated = ct.electron(_func, **decorator_kwargs)
        elif wflow_engine == "jobflow":
            from jobflow import job as jf_job

            decorated = jf_job(_func, **decorator_kwargs)
        elif wflow_engine == "parsl":
            from parsl import python_app

            decorated = python_app(_func, **decorator_kwargs)
        else:
            decorated = _func

        return decorated(*f_args, **f_kwargs)

    if _func is None:

        def decorator(_f):
            return job(_f, **kwargs)

        return decorator

    return _inner


def flow(
    _func: callable | None = None,
    engine: Literal["covalent", "parsl", "jobflow", "local"] | None = None,
    **kwargs,
) -> Flow:  # sourcery skip
    """
    Decorator for workflows, which consist of at least one compute job. This is a @flow decorator.

    flow = ct.lattice/ct.dispatch(ct.lattice) [Covalent]. For Parsl and Jobflow, the decorator returns the
    original function, unchanged.

    The wrapped function gets a new kwarg, `decorator_kwargs`, that can
    be used to modify the workflow engine decorator keyword arguments.
    The wrapped function can also be accessed without the decorator
    by calling the `.__wrapped__` attribute.

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    engine
        The workflow engine. This defaults to `SETTINGS.WORKFLOW_ENGINE` if `None`.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    Flow
        The @flow-decorated function.
    """

    @functools.wraps(_func)
    def _inner(
        *f_args,
        decorator_kwargs: dict | None = None,
        dispatch_kwargs: dict | None = None,
        **f_kwargs,
    ) -> Any:
        """
        The @flow-decorated function.

        Parameters
        ----------
        *f_args
            Positional arguments to the function, if any.
        decorator_kwargs
            Keyword arguments to pass to the workflow engine decorator.
        dispatch_kwargs
            Keyword arguments to pass to ct.dispatch() if Covalent is used.
        **f_kwargs
            Keyword arguments to the function, if any.

        Returns
        -------
        Any
            The output of the @flow-decorated function.
        """
        from quacc import SETTINGS

        if not decorator_kwargs:
            decorator_kwargs = kwargs

        wflow_engine = engine or SETTINGS.WORKFLOW_ENGINE

        dispatch_kwargs = dispatch_kwargs or {}
        if dispatch_kwargs and wflow_engine != "covalent":
            raise ValueError("The `dispatch_kwargs` argument only works with Covalent.")

        if wflow_engine == "covalent":
            import covalent as ct

            try:
                return ct.dispatch(ct.lattice(_func, **decorator_kwargs))(
                    *f_args, **f_kwargs
                )
            except (AttributeError, TypeError):
                return ct.lattice(_func, **dispatch_kwargs, **decorator_kwargs)(
                    *f_args, **f_kwargs
                )
        else:
            decorated = _func

        return decorated(*f_args, **f_kwargs)

    if _func is None:

        def decorator(_f):
            return flow(_f, **kwargs)

        return decorator

    return _inner


def subflow(
    _func: callable | None = None,
    engine: Literal["covalent", "parsl", "jobflow", "local"] | None = None,
    **kwargs,
) -> Subflow:  # sourcery skip
    """
    Decorator for (dynamic) sub-workflows. This is a @subflow decorator.

    subflow = ct.electron(ct.lattice) [Covalent] = join_app [Parsl].
    For Jobflow, the decorator returns the original (unwrapped) function.

    The wrapped function gets a new kwarg, `decorator_kwargs`, that can
    be used to modify the workflow engine decorator keyword arguments.
    The wrapped function can also be accessed without the decorator
    by calling the `.__wrapped__` attribute.

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    engine
        The workflow engine. This defaults to `SETTINGS.WORKFLOW_ENGINE` if `None`.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    callable
        The decorated function.
    """

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        """
        The @subflow-decorated function.

        Parameters
        ----------
        *f_args
            Positional arguments to the function, if any.
        decorator_kwargs
            Keyword arguments to pass to the workflow engine decorator.
        **f_kwargs
            Keyword arguments to the function, if any.

        Returns
        -------
        Any
            The output of the @subflow-decorated function.
        """
        from quacc import SETTINGS

        if not decorator_kwargs:
            decorator_kwargs = kwargs

        wflow_engine = engine or SETTINGS.WORKFLOW_ENGINE
        if wflow_engine == "covalent":
            import covalent as ct

            decorated = ct.electron(ct.lattice(_func), **decorator_kwargs)
        elif wflow_engine == "parsl":
            from parsl import join_app

            decorated = join_app(_func, **decorator_kwargs)
        else:
            decorated = _func

        return decorated(*f_args, **f_kwargs)

    if _func is None:

        def decorator(_f):
            return subflow(_f, **kwargs)

        return decorator

    return _inner


def fetch_atoms(atoms: Atoms | dict) -> Atoms:
    """
    Returns an Atoms object from a typical quacc recipe input, which can
    either be an `Atoms` object or a dictionary with an entry `{"atoms": Atoms}`.
    It may seem a bit silly to do this, but there is a purpose. If you want to
    write a workflow where the output of one recipe is passed to the input of
    another recipe, you can always do output["atoms"] to fetch the output Atoms
    object to pass to the input to the second function. However, this process
    will often be its own compute step in workflow management tools because they
    need to resolve the output in order to query it. Depending on the workflow manager,
    this can be a waste of compute resources, so it's oftentimes better to do this
    parsing inside the compute task itself, which is why passing in the raw dictionary
    can be preferred.

    Parameters
    ----------
    atoms
        Atoms object or dictionary with an entry {"atoms": Atoms}

    Returns
    -------
    Atoms
        Atoms object
    """
    try:
        return atoms["atoms"]
    except Exception:
        return atoms
