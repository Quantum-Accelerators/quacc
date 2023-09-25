import pytest

try:
    from quacc.utils.wflows import _make_dask_cluster, make_prefect_runner

    from dask_jobqueue import PBSCluster, SLURMCluster  # isort: skip
    from prefect_dask.task_runners import DaskTaskRunner  # isort: skip

    dask_prefect = True
except ImportError:
    dask_prefect = None

pytestmark = pytest.mark.skipif(
    dask_prefect is None,
    reason="Dask and Prefect dependencies must be installed.",
)


def test_make_prefect_runner():
    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    runner = make_prefect_runner(cluster_kwargs, temporary=True)
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_prefect_runner(
        cluster_kwargs, adapt_kwargs=adapt_kwargs, temporary=True
    )
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
    runner = make_prefect_runner(
        cluster_kwargs, adapt_kwargs=adapt_kwargs, temporary=True
    )
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == SLURMCluster
    assert runner.cluster_kwargs == cluster_kwargs
    assert runner.adapt_kwargs == adapt_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    runner = make_prefect_runner(
        cluster_kwargs, cluster_class=PBSCluster, temporary=True
    )
    assert isinstance(runner, DaskTaskRunner)
    assert runner.cluster_class == PBSCluster
    assert runner.cluster_kwargs == cluster_kwargs

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    adapt_kwargs = {"minimum": 1, "maximum": 2}
    runner = make_prefect_runner(cluster_kwargs, adapt_kwargs=adapt_kwargs)


def test_make_dask_cluster():
    from dask_jobqueue import SLURMCluster

    cluster_kwargs = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
    }
    cluster = _make_dask_cluster(SLURMCluster, cluster_kwargs, verbose=True)
    assert isinstance(cluster, SLURMCluster)
