import pytest

try:
    from dask_jobqueue import PBSCluster, SLURMCluster
    from prefect_dask.task_runners import DaskTaskRunner

    from quacc.util.dask import _make_cluster, make_runner

    dask_prefect = True
except ImportError:
    dask_prefect = None


@pytest.mark.skipif(
    dask_prefect is None,
    reason="Dask and Prefect dependencies must be installed.",
)
def test_make_runner():
    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    runner = make_runner(cluster_kwargs, temporary=True)
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_runner(cluster_kwargs, adapt_kwargs=adapt_kwargs, temporary=True)
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs
    assert runner.adapt_kwargs == adapt_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_runner(cluster_kwargs, adapt_kwargs=adapt_kwargs, temporary=True)
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs
    assert runner.adapt_kwargs == adapt_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    runner = make_runner(cluster_kwargs, cluster_class=PBSCluster, temporary=True)
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == PBSCluster
    assert runner.cluster_kwargs == cluster_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    runner = make_runner(cluster_kwargs)


@pytest.mark.skipif(
    dask_prefect is None,
    reason="Dask and Prefect dependencies must be installed.",
)
def test_make_cluster():
    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    cluster = _make_cluster(SLURMCluster, cluster_kwargs, verbose=True)
    assert isinstance(cluster, SLURMCluster)
