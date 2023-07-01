"""
Utility functions for workflow engines
"""
from __future__ import annotations

import asyncio

from dask_jobqueue import SLURMCluster
from dask_jobqueue.core import Job


def make_dask_cluster(
    cluster_params: dict, cluster: Job = SLURMCluster, num_jobs: int = 1
) -> Job:
    """
    Make a Dask cluster for use with Prefect workflows.

    Parameters
    ----------
    cluster_params
        Keyword arguments to pass to `cluster`.
    cluster
        The Dask cluster to use. Defaults to `SLURMCluster`.
    n_jobs
        Number of Slurm jobs to run on the Dask cluster.
    """

    async def make_cluster(num_jobs: int, cluster: Job, cluster_params: dict) -> Job:
        custom_cluster = await cluster(**cluster_params)
        custom_cluster.scale(num_jobs)
        return custom_cluster

    return asyncio.run(make_cluster(num_jobs, cluster, cluster_params))
