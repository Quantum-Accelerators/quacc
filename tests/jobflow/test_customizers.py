from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")
from pathlib import Path

from quacc import get_settings, job, redecorate, strip_decorator


def test_strip_decorators():
    @job
    def add(a, b):
        return a + b

    stripped_add = strip_decorator(add)
    assert stripped_add(1, 2) == 3


def test_change_settings_redecorate_job(tmp_path_factory, jobflow_output):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")

    @job
    def write_file_job(name="job.txt"):
        path = Path(get_settings().SCRATCH_DIR, name)
        with open(path, "w") as f:
            f.write("test file")
        return path

    write_file_job = redecorate(
        write_file_job, job(settings_swap={"SCRATCH_DIR": tmp_dir1})
    )
    job = write_file_job()
    responses = jf.run_locally(job, ensure_success=True, create_folders=False)
    assert Path(jobflow_output(responses, job)) == tmp_dir1 / "job.txt"
    assert Path(tmp_dir1 / "job.txt").exists()
