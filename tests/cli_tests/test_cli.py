from subprocess import PIPE, Popen

import covalent as ct


def test_config():
    process = Popen(["quacc", "config"], stdout=PIPE)
    _, stderr = process.communicate()
    assert stderr is None
    ct_config = ct.get_config()
    for executor in ct_config["executors"]:
        if "create_unique_workdir" in ct_config["executors"][executor]:
            assert ct_config["executors"][executor]["create_unique_workdir"] is True
    if "slurm" in ct_config["executors"]:
        assert ct_config["executors"]["slurm"].get("use_srun", True) is False


def test_help():
    process = Popen("quacc", stdout=PIPE)
    stdout, stderr = process.communicate()
    assert stderr is None
    assert "Show this message and exit" in stdout.decode("utf-8")
