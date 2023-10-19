import os

from quacc.utils.files import make_unique_dir


def test_make_unique_dir(tmpdir):
    tmpdir.chdir()

    jobdir = make_unique_dir()
    assert os.path.exists(jobdir)

    jobdir = make_unique_dir(base_path="tmp_dir")
    assert os.path.exists("tmp_dir")
    assert "tmp_dir" in jobdir
    assert os.path.exists(jobdir)
