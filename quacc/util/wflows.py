"""
Utility functions for workflow engines
"""
from __future__ import annotations

import asyncio

from dask_jobqueue import SLURMCluster
from dask_jobqueue.core import Job
from prefect_dask.task_runners import DaskTaskRunner


def make_dask_cluster(
    cluster_params: dict, num_jobs: int = 1, cluster: Job = SLURMCluster
) -> DaskTaskRunner:
    """
    Spin up a Dask cluster for use with Prefect workflows.

    Parameters
    ----------
    cluster_params
        Keyword arguments to pass to `cluster`.
    n_jobs
        Number of Slurm jobs to run on the Dask cluster.
    cluster
        The Dask cluster to use. Defaults to `SLURMCluster`.

    Returns
    -------
    DaskTaskRunner
        A DaskTaskRunner object for use with Prefect workflows.
    """

    async def make_cluster(
        cluster_params: dict,
        num_jobs: int,
        cluster: Job,
    ) -> Job:
        """
        Make a Dask cluster.

        Parameters
        ----------
        cluster_params
            Keyword arguments to pass to `cluster`.
        num_jobs
            Number of Slurm jobs to run on the Dask cluster.
        cluster
            The Dask cluster to use. Defaults to `SLURMCluster`.
        """
        custom_cluster = await cluster(**cluster_params)
        custom_cluster.scale(num_jobs)
        return custom_cluster

    asyncio.run(make_cluster(cluster_params, num_jobs, cluster))

    return DaskTaskRunner(cluster.scheduler_address)
