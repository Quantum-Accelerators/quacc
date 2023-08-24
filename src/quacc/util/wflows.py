from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from covalent import electron as ct_electron
    from covalent import lattice as ct_lattice
    from dask_jobqueue.core import DaskJobqueueJob
    from jobflow import Job as JobflowJob
    from parsl.app.python import PythonApp


def job(
    _func: callable | None = None, **kwargs
) -> callable | ct_electron | JobflowJob | PythonApp:
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
        The decorated function. The decorated function will have an attribute `original_func`
        which is the undecorated function.
    """

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.electron(_func, **kwargs)
    elif wflow_engine == "jobflow":
        from jobflow import job as jf_job

        decorated = jf_job(_func, **kwargs)
    elif wflow_engine == "parsl":
        from parsl import python_app

        decorated = python_app(_func, **kwargs)
    elif not wflow_engine:
        decorated = _func
    else:
        msg = f"Unknown workflow engine: {wflow_engine}"
        raise ValueError(msg)

    decorated.original_func = _func

    return decorated


def flow(_func: callable | None = None, **kwargs) -> callable | ct_lattice:
    """
    Decorator for workflows, which consist of at least one compute job. This is a @flow decorator.

    @flow = @ct.lattice [Covalent]. For Parsl and Jobflow, the decorator returns the
    undecorated function.

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

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.lattice(_func, **kwargs)
    elif wflow_engine in {"jobflow", "parsl"} or not wflow_engine:
        decorated = _func
    else:
        msg = f"Unknown workflow engine: {wflow_engine}"
        raise ValueError(msg)

    return decorated


def subflow(_func: callable | None = None, **kwargs) -> callable | ct_electron:
    """
    Decorator for (dynamic) sub-workflows. This is a @subflow decorator.

    @subflow = @ct.electron(@ct.lattice) [Covalent] = @join_app [Parsl].
    For Jobflow, the decorator returns the undecorated function.

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

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.electron(ct.lattice(_func), **kwargs)
    elif wflow_engine == "parsl":
        from parsl import join_app

        decorated = join_app(_func, **kwargs)
    elif wflow_engine == "jobflow" or not wflow_engine:
        decorated = _func
    else:
        msg = f"Unknown workflow engine: {wflow_engine}"
        raise ValueError(msg)

    return decorated
