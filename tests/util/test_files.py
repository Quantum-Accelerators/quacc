import os
from shutil import rmtree

from quacc.util.files import make_unique_dir


def teardown_function():
    for f in os.listdir(os.getcwd()):
        if "quacc_" in f or "tmp_dir" in f:
            rmtree(f)


def test_make_unique_dir():
    jobdir = make_unique_dir()
    assert os.path.exists(jobdir)

    jobdir = make_unique_dir(base_path="tmp_dir")
    assert os.path.exists("tmp_dir")
    assert "tmp_dir" in jobdir
    assert os.path.exists(jobdir)
