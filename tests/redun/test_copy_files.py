from __future__ import annotations

import pytest

redun = pytest.importorskip("redun")

from pathlib import Path

from quacc import flow, job
from quacc.wflow_tools.job_argument import Copy


@pytest.fixture
def scheduler():
    return redun.Scheduler()


def test_copy_files(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    @job
    def create_file(name: str, copy: Copy | None = None):
        output_dir = tmp_path / name
        output_dir.mkdir(parents=True, exist_ok=True)
        Path(output_dir / name).touch()

        extra = None
        if copy is not None:
            extra = copy.do_copy(output_dir)

        # Note: Important to return the copy.do_copy task for it to be run!
        return {"dir_name": output_dir, "extra": extra}

    @flow
    def create_files():
        job1 = create_file("job1")
        job2 = create_file("job2", copy=Copy({job1["dir_name"]: "job1*"}))
        return [job1, job2]

    scheduler.run(create_files())

    # Individual job folders/files should exist, and the job1 file should be
    # copied over to the job2 folder.
    assert Path(tmp_path / "job1/job1").exists()
    assert Path(tmp_path / "job2/job2").exists()
    assert Path(tmp_path / "job2/job1").exists()
