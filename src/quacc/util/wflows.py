from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

if TYPE_CHECKING:
    from covalent import electron as ct_electron
    from covalent import lattice as ct_lattice
    from jobflow import Job as JobflowJob
    from parsl.app.python import PythonApp
    from prefect import Flow as PrefectFlow
    from prefect import Task as PrefectTask


try:
    from dask_jobqueue import SLURMCluster
    from prefect_dask.task_runners import DaskTaskRunner

    if TYPE_CHECKING:
        from dask_jobqueue.core import DaskJobqueueJob

    prefect_deps = True

except ImportError:
    prefect_deps = False


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
        The decorated function.
    """

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        return ct.electron(_func, **kwargs)
    if wflow_engine == "jobflow":
        from jobflow import job as jf_job

        return jf_job(_func, **kwargs)
    if wflow_engine == "parsl":
        from parsl import python_app

        return python_app(_func, **kwargs)
    if wflow_engine == "prefect":
        from prefect import task

        return task(_func, **kwargs)

    return _func


def flow(
    _func: callable | None = None, **kwargs
) -> callable | ct_lattice | PrefectFlow:
    """
    Decorator for workflows, which consist of at least one compute job. This is a @flow decorator.

    @flow = @ct.lattice [Covalent] = @flow [Prefect]. For Prefect, the decorator returns the
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

        return ct.lattice(_func, **kwargs)
    if wflow_engine == "jobflow":
        raise NotImplementedError(
            "Jobflow is not compatible with the use of a @flow decorator. Instead, you should use the `Flow()` object in Jobflow to stitch together individual compute jobs."
        )
    if wflow_engine == "parsl":
        return _func
    if wflow_engine == "prefect":
        from prefect import flow as prefect_flow

        return prefect_flow(_func, **kwargs)

    return _func


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

        return ct.electron(ct.lattice(_func), **kwargs)
    if wflow_engine == "jobflow":
        raise NotImplementedError(
            "Jobflow is not compatible with the use of a @subflow decorator. Instead, you should use the `Response` object in Jobflow to create a dynamic workflow."
        )
    if wflow_engine == "parsl":
        from parsl import join_app

        return join_app(_func, **kwargs)
    if wflow_engine == "prefect":
        from prefect import flow as prefect_flow

        return prefect_flow(_func, **kwargs)

    return _func


@requires(prefect_deps, "Need quacc[prefect] dependencies")
def make_runner(
    cluster_kwargs: dict,
    cluster_class: callable = None,
    adapt_kwargs: dict[str, int | None] | None = None,
    client_kwargs: dict = None,
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
    cluster = _make_cluster(cluster_class, cluster_kwargs)

    # Set up adaptive scaling
    if adapt_kwargs and (adapt_kwargs["minimum"] or adapt_kwargs["maximum"]):
        cluster.adapt(minimum=adapt_kwargs["minimum"], maximum=adapt_kwargs["maximum"])

    # Return the DaskTaskRunner with the cluster address
    return DaskTaskRunner(address=cluster.scheduler_address)


@requires(prefect_deps, "Need quacc[prefect] dependencies")
def _make_cluster(
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
