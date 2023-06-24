import os

import covalent as ct


def test_config():
    os.system("quacc config")
    ct_config = ct.get_config()
    for executor in ct_config["executors"]:
        if "create_unique_workdir" in ct_config["executors"][executor]:
            assert ct_config["executors"][executor]["create_unique_workdir"] is True
    if "slurm" in ct_config["executors"]:
        assert ct_config["executors"]["slurm"].get("use_srun", True) is False


def test_help():
    os.system("quacc help")
