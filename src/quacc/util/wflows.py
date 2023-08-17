from __future__ import annotations

from monty.dev import requires

from quacc import SETTINGS

try:
    from dask_jobqueue import SLURMCluster
    from prefect_dask.task_runners import DaskTaskRunner

    if TYPE_CHECKING:
        from dask_jobqueue.core import Job

    prefect_deps = True

except ImportError:
    prefect_deps = False


def job(_func: callable | None = None, **kwargs):
    """
    Decorator for workflow jobs
    """

    wflow_manager = (
        SETTINGS.WORKFLOW_MANAGER.lower() if SETTINGS.WORKFLOW_MANAGER else None
    )

    def wrapper(*func_args, **func_kwargs):
        if not wflow_manager:
            return _func(*func_args, **func_kwargs)
        elif wflow_manager == "covalent":
            import covalent as ct

            return ct.electron(_func, **kwargs)(*func_args, **func_kwargs)
        elif wflow_manager == "parsl":
            from parsl import python_app

            return python_app(_func, **kwargs)(*func_args, **func_kwargs)
        elif wflow_manager == "prefect":
            from prefect import task

            return task(_func, **kwargs)(*func_args, **func_kwargs)
        elif wflow_manager == "jobflow":
            from jobflow import job as jf_job

            return jf_job(_func, **kwargs)(*func_args, **func_kwargs)
        else:
            raise ValueError(f"Unknown workflow manager {wflow_manager}.")

    return wrapper


def flow(_func: callable | None = None, **kwargs):
    """
    Decorator for workflow flows
    """

    wflow_manager = (
        SETTINGS.WORKFLOW_MANAGER.lower() if SETTINGS.WORKFLOW_MANAGER else None
    )

    def wrapper(*func_args, **func_kwargs):
        if not wflow_manager or wflow_manager == "parsl":
            return _func(*func_args, **func_kwargs)
        elif wflow_manager == "covalent":
            import covalent as ct

            return ct.lattice(_func, **kwargs)(*func_args, **func_kwargs)
        elif wflow_manager == "prefect":
            from prefect import flow

            return flow(_func, **kwargs)(*func_args, **func_kwargs)
        else:
            raise ValueError(f"Unknown workflow manager {wflow_manager}.")

    return wrapper


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
def _make_cluster(cluster_class: callable, cluster_kwargs: dict, verbose=True) -> Job:
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
