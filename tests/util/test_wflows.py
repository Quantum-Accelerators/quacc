import pytest

try:
    from quacc.util.wflows import make_dask_cluster

    dask = True
except ImportError:
    dask = None


@pytest.mark.skipif(
    dask is None,
    reason="dask-jobqueue must be installed.",
)
def test_make_cluster():
    cluster_params = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
        "account": "MyAccount",
    }
    make_dask_cluster(cluster_params)
