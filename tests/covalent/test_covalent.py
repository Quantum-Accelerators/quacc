import pytest
import quacc

try:
    import covalent as ct
except ImportError:
    ct = None


@pytest.mark.skipif(
    ct is None,
    reason="covalent must be installed. Try pip install covalent",
)
def test_covalent_config():
    ct_config = ct.get_config()
    for executor in ct_config["executors"]:
        assert (
            ct_config["executors"][executor].get("create_unique_workdir", "false")
            == "true"
        )
    if ct_config["executors"].get("slurm", None):
        assert ct_config["executors"]["slurm"].get("use_srun", "true") == "false"
