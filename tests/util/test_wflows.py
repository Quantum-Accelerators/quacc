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
def test_make_cluster(capsys):
    # NOTE: I don't know how to either turn off the stderr
    # output from dask-jobqueue or mock it.
    cluster_params = {
        "cores": 1,
        "memory": "1GB",
        "processes": 1,
        "account": "MyAccount",
    }
    with capsys.disabled():
        make_dask_cluster(cluster_params)
