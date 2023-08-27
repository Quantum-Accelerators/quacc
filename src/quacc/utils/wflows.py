from __future__ import annotations

import functools
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms


def job(_func: callable | None = None, **kwargs) -> callable:
    """
    Decorator for individual compute jobs. This is a @job decorator.

    @job = @ct.electron [Covalent] = @python_app [Parsl] = @job [Jobflow]

    Parameters
    ----------
    _func
        The function to decorate.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    callable
        The decorated function. The decorated function will have an attribute `__wrapped__`
        which is the __wrapped__ function.
    """

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        from quacc import SETTINGS

        if decorator_kwargs is None:
            decorator_kwargs = kwargs or {}

        wflow_engine = (
            SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
        )
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

    return _inner


def flow(_func: callable | None = None, **kwargs) -> callable:
    """
    Decorator for workflows, which consist of at least one compute job. This is a @flow decorator.

    @flow = @ct.lattice [Covalent]. For Parsl and Jobflow, the decorator returns the
    __wrapped__ function.

    Parameters
    ----------
    _func
        The function to decorate.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    callable
        The decorated function.
    """

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        from quacc import SETTINGS

        if decorator_kwargs is None:
            decorator_kwargs = kwargs or {}

        wflow_engine = (
            SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
        )
        if wflow_engine == "covalent":
            import covalent as ct

            decorated = ct.lattice(_func, **decorator_kwargs)
        else:
            decorated = _func

        return decorated(*f_args, **f_kwargs)

    return _inner


def subflow(_func: callable | None = None, **kwargs) -> callable:
    """
    Decorator for (dynamic) sub-workflows. This is a @subflow decorator.

    @subflow = @ct.electron(@ct.lattice) [Covalent] = @join_app [Parsl].
    For Jobflow, the decorator returns the __wrapped__ function.

    Parameters
    ----------
    _func
        The function to decorate.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    callable
        The decorated function.
    """

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        from quacc import SETTINGS

        if decorator_kwargs is None:
            decorator_kwargs = kwargs or {}

        wflow_engine = (
            SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
        )
        if wflow_engine == "covalent":
            import covalent as ct

            decorated = ct.electron(ct.lattice(_func), **decorator_kwargs)
        elif wflow_engine == "parsl":
            from parsl import join_app

            decorated = join_app(_func, **decorator_kwargs)
        else:
            decorated = _func

        return decorated(*f_args, *f_kwargs)

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
    need to resolve the output in order to query it. For instance, in Parsl you'd
    have to call `.result()` and start a block.

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
