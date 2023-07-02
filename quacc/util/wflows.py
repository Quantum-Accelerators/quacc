"""Utility functions for workflow engines"""
from __future__ import annotations

import asyncio

from monty.dev import requires

try:
    from dask_jobqueue import SLURMCluster
    from dask_jobqueue.core import Job
    from prefect_dask.task_runners import DaskTaskRunner

except:
    prefect = None


@requires(prefect, "Install quacc[prefect] extras")
def launch_runner(
    cluster_kwargs: dict, cluster_class: callable = SLURMCluster, verbose: bool = False
) -> DaskTaskRunner:
    """
    Make a DaskTaskRunner for use with Prefect workflows. This function will immediately
    submit a Slurm job upon being called and will wait for work to be run.

    Parameters
    ----------
    cluster_kwargs
        Keyword arguments to pass to `cluster`.
    cluster_class
        The Dask cluster class to use. Defaults to `dask_jobqueue.SLURMCluster`.
    verbose
        Whether to print out the job script or not (useful for debugging purposes).

    Returns
    -------
    DaskTaskRunner
        A DaskTaskRunner object for use with Prefect workflows.
    """

    async def make_cluster(cluster_class: callable, cluster_kwargs: dict) -> Job:
        cluster = await cluster_class(**cluster_kwargs)
        if verbose:
            print(cluster.job_script())
        return cluster

    cluster = asyncio.run(make_cluster(cluster_class, cluster_kwargs))

    return DaskTaskRunner(cluster.scheduler_address)


@requires(prefect, "Install quacc[prefect] extras")
def make_runner(
    cluster_kwargs: dict,
    cluster_class: callable = SLURMCluster,
) -> DaskTaskRunner:
    """
    Make a DaskTaskRunner for use with Prefect workflows. This DaskTaskRunner
    will only submit a Slurm job once the `Flow` begins.

    Parameters
    ----------
    cluster_kwargs
        Keyword arguments to pass to `cluster_class`.
    cluster_class
        The Dask cluster class to use. Defaults to `dask_jobqueue.SLURMCluster`.

    Returns
    -------
    DaskTaskRunner
        A DaskTaskRunner object for use with Prefect workflows.
    """
    return DaskTaskRunner(cluster_class=cluster_class, cluster_kwargs=cluster_kwargs)
