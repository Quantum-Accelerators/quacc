import pytest

from quacc.wflow_tools.prefect import _make_dask_cluster, make_prefect_runner

pytest.importorskip("dask_jobqueue")
pytest.importorskip("prefect_dask")

from dask_jobqueue import PBSCluster, SLURMCluster
from prefect_dask.task_runners import DaskTaskRunner


def test_make_prefect_runner():
    cluster_kwargs = {"cores": 1, "memory": "1GB", "processes": 1}
    runner = make_prefect_runner(cluster_kwargs, temporary=True)
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs

    cluster_kwargs = {"cores": 1, "memory": "1GB", "processes": 1}
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_prefect_runner(
        cluster_kwargs, adapt_kwargs=adapt_kwargs, temporary=True
    )
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs
    assert runner.adapt_kwargs == adapt_kwargs

    cluster_kwargs = {"cores": 1, "memory": "1GB", "processes": 1}
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_prefect_runner(
        cluster_kwargs, adapt_kwargs=adapt_kwargs, temporary=True
    )
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs
    assert runner.adapt_kwargs == adapt_kwargs

    cluster_kwargs = {"cores": 1, "memory": "1GB", "processes": 1}
    runner = make_prefect_runner(
        cluster_kwargs, cluster_class=PBSCluster, temporary=True
    )
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == PBSCluster
    assert runner.cluster_kwargs == cluster_kwargs

    cluster_kwargs = {"cores": 1, "memory": "1GB", "processes": 1}
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_prefect_runner(cluster_kwargs, adapt_kwargs=adapt_kwargs)


def test_make_dask_cluster():
    cluster_kwargs = {"cores": 1, "memory": "1GB", "processes": 1}
    cluster = _make_dask_cluster(SLURMCluster, cluster_kwargs, verbose=True)
    assert isinstance(cluster, SLURMCluster)
