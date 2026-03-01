from __future__ import annotations

import pytest

prefect = pytest.importorskip("prefect")

from pathlib import Path

from quacc import flow, job
from quacc.wflow_tools.job_argument import Copy


def test_copy_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def create_file(name: str, copy: Copy | None = None):
        output_dir = tmp_path / name
        output_dir.mkdir(parents=True, exist_ok=True)
        Path(output_dir / name).touch()

        if copy is not None:
            copy.do_copy(output_dir)

        return {"dir_name": output_dir}

    @flow
    def create_files():
        job1 = create_file("job1")
        job2 = create_file("job2", copy=Copy({job1["dir_name"]: "job1*"}))
        return [job1, job2]

    create_files()

    # Individual job folders/files should exist, and the job1 file should be
    # copied over to the job2 folder.
    assert Path(tmp_path / "job1/job1").exists()
    assert Path(tmp_path / "job2/job2").exists()
    assert Path(tmp_path / "job2/job1").exists()
