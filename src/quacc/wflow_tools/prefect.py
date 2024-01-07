"""Utilities for Prefect."""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

try:
    from prefect_dask.task_runners import DaskTaskRunner

    prefect_deps = True
except ImportError:
    prefect_deps = False
try:
    from dask_jobqueue import SLURMCluster

    dask_deps = True
except ImportError:
    dask_deps = False

if TYPE_CHECKING:
    from typing import Any, Callable

    from dask_jobqueue.core import Job as DaskJob


@requires(prefect_deps and dask_deps, "Need quacc[prefect] dependencies")
def make_prefect_runner(
    cluster_kwargs: dict[str, Any],
    cluster_class: Callable | None = None,
    adapt_kwargs: dict[str, int | None] | None = None,
    client_kwargs: dict[str, Any] | None = None,
    temporary: bool = False,
) -> DaskTaskRunner:
    """
    Make a `DaskTaskRunner` for use with Prefect workflows.

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
    cluster = _make_dask_cluster(cluster_class, cluster_kwargs)

    # Set up adaptive scaling
    if adapt_kwargs and (adapt_kwargs["minimum"] or adapt_kwargs["maximum"]):
        cluster.adapt(minimum=adapt_kwargs["minimum"], maximum=adapt_kwargs["maximum"])

    # Return the DaskTaskRunner with the cluster address
    return DaskTaskRunner(address=cluster.scheduler_address)


@requires(dask_deps, "Need quacc[prefect] dependencies")
def _make_dask_cluster(
    cluster_class: Callable, cluster_kwargs: dict[str, Any], verbose=False
) -> DaskJob:
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
        print(cluster.job_script())  # noqa: T201
    return cluster
