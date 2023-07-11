from __future__ import annotations

from typing import Literal

try:
    from dask_jobqueue import SLURMCluster
    from dask_jobqueue.core import Job
    from prefect_dask.task_runners import DaskTaskRunner

except ImportError:
    raise ImportError("Install quacc[prefect]")


def make_runner(
    cluster_kwargs: dict,
    cluster_class: callable = SLURMCluster,
    adapt_kwargs: dict[str, int | None] | None = None,
    client_kwargs: dict = None,
    temporary: bool = True,
) -> DaskTaskRunner:
    """
    Make a DaskTaskRunner for use with Prefect workflows. This DaskTaskRunner
    will only submit a Slurm job once the `Flow` begins. It will terminate
    once the `Flow` is finished.

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


def _make_cluster(cluster_class: callable, cluster_kwargs: dict, verbose=False) -> Job:
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
        print(cluster.job_script())
    return cluster
