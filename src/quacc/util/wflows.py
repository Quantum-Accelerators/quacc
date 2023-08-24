from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from covalent import electron as ct_electron
    from covalent import lattice as ct_lattice
    from dask_jobqueue.core import DaskJobqueueJob
    from jobflow import Job as JobflowJob
    from parsl.app.python import PythonApp
    from prefect import Flow as PrefectFlow
    from prefect import Task as PrefectTask
    from prefect_dask.task_runners import DaskTaskRunner


def job(
    _func: callable | None = None, **kwargs
) -> callable | ct_electron | JobflowJob | PythonApp | PrefectTask:
    """
    Decorator for individual compute jobs. This is a @job decorator.

    @job = @ct.electron [Covalent] = @python_app [Parsl] = @job [Jobflow] = @task [Prefect]

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
    elif wflow_engine == "prefect":
        from prefect import task

        decorated = task(_func, **kwargs)
    elif not wflow_engine:
        decorated = _func
    else:
        msg = f"Unknown workflow engine: {wflow_engine}"
        raise ValueError(msg)

    decorated.original_func = _func

    return decorated


def flow(
    _func: callable | None = None, **kwargs
) -> callable | ct_lattice | PrefectFlow:
    """
    Decorator for workflows, which consist of at least one compute job. This is a @flow decorator.

    @flow = @ct.lattice [Covalent] = @flow [Prefect]. For Parsl, the decorator returns the
    undecorated function. This decorator is not compatible with jobflow.

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
    elif wflow_engine == "prefect":
        from prefect import flow as prefect_flow

        decorated = prefect_flow(_func, **kwargs)
    elif wflow_engine in {"jobflow", "parsl"} or not wflow_engine:
        decorated = _func
    else:
        msg = f"Unknown workflow engine: {wflow_engine}"
        raise ValueError(msg)

    return decorated


def subflow(
    _func: callable | None = None, **kwargs
) -> callable | ct_electron | PrefectFlow:
    """
    Decorator for (dynamic) sub-workflows. This is a @subflow decorator.

    @subflow = @ct.electron(@ct.lattice) [Covalent] = @join_app [Parsl] = @flow [Prefect].
    This decorator is not compatible with jobflow.

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
    elif wflow_engine == "prefect":
        from prefect import flow as prefect_flow

        decorated = prefect_flow(_func, **kwargs)
    elif wflow_engine == "jobflow" or not wflow_engine:
        decorated = _func
    else:
        msg = f"Unknown workflow engine: {wflow_engine}"
        raise ValueError(msg)

    return decorated


def make_dask_runner(
    cluster_kwargs: dict,
    cluster_class: callable | None = None,
    adapt_kwargs: dict[str, int | None] | None = None,
    client_kwargs: dict | None = None,
    temporary: bool = False,
) -> DaskTaskRunner:
    """
    Make a DaskTaskRunner for use with Prefect workflows.

    Parameters
    ----------
    cluster_kwargs
        Keyword arguments to pass to `cluster_class`.
    cluster_class
        The Dask cluster class to use. Defaults to `dask_jobqueue.SLURMCluster`.
    adapt_kwargs
        Keyword arguments to pass to `cluster.adapt` of the form `{"minimum": int, "maximum": int}`.
        If `None`, no adaptive scaling will be done.
    client_kwargs
        Keyword arguments to pass to `dask.distributed.Client`.
    temporary
        Whether to use a temporary cluster. If `True`, the cluster will be
        terminated once the `Flow` is finished. If `False`, the cluster will
        run until the walltime is reached and can run multiple `Flow`s.

    Returns
    -------
    DaskTaskRunner
        A DaskTaskRunner object for use with Prefect workflows.
    """
    from dask_jobqueue import SLURMCluster
    from prefect_dask.task_runners import DaskTaskRunner

    if cluster_class is None:
        cluster_class = SLURMCluster

    # Make the one-time-use DaskTaskRunner
    if temporary:
        return DaskTaskRunner(
            cluster_class=cluster_class,
            cluster_kwargs=cluster_kwargs,
            adapt_kwargs=adapt_kwargs,
            client_kwargs=client_kwargs,
        )

    # Make the Dask cluster
    cluster = _make_dask_cluster(cluster_class, cluster_kwargs)

    # Set up adaptive scaling
    if adapt_kwargs and (adapt_kwargs["minimum"] or adapt_kwargs["maximum"]):
        cluster.adapt(minimum=adapt_kwargs["minimum"], maximum=adapt_kwargs["maximum"])

    # Return the DaskTaskRunner with the cluster address
    return DaskTaskRunner(address=cluster.scheduler_address)


def _make_dask_cluster(
    cluster_class: callable, cluster_kwargs: dict, verbose=True
) -> DaskJobqueueJob:
    """
    Make a Dask cluster for use with Prefect workflows.

    Parameters
    ----------
    cluster_class
        The Dask cluster class to use. Defaults to `dask_jobqueue.SLURMCluster`.
    cluster_kwargs
        Keyword arguments to pass to `cluster_class`.
    verbose
        Whether to print the job script to stdout.
    """

    cluster = cluster_class(**cluster_kwargs)
    if verbose:
        print(
            f"Workers are submitted with the following job script:\n{cluster.job_script()}"
        )
        print(f"Scheduler is running at {cluster.scheduler.address}")
        print(f"Dashboard is located at {cluster.dashboard_link}")

    return cluster
