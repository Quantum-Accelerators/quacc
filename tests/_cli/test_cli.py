import covalent as ct
from click.testing import CliRunner

from quacc._cli.cli import cli


def test_help():
    response = CliRunner().invoke(cli).output
    assert "Show this message and exit" in response


def test_config():
    response = CliRunner().invoke(cli, "config").output
    assert "Covalent configuration: Complete" in response
    ct_config = ct.get_config()
    for executor in ct_config["executors"]:
        if "create_unique_workdir" in ct_config["executors"][executor]:
            assert ct_config["executors"][executor]["create_unique_workdir"] is True
    if "slurm" in ct_config["executors"]:
        assert ct_config["executors"]["slurm"].get("use_srun", True) is False
