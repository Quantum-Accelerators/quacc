from __future__ import annotations

import pytest

redun = pytest.importorskip("redun")

from pathlib import Path

from quacc import flow, job
from quacc.utils.files import copy_decompress_files


@pytest.fixture
def scheduler():
    return redun.Scheduler()


def test_copy_files(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)

    @job
    def create_file(name: str, copy: list[dict] | None = None):
        output_dir = tmp_path / name
        output_dir.mkdir(parents=True, exist_ok=True)
        Path(output_dir / name).touch()

        if copy is not None:
            for spec in copy:
                copy_decompress_files(spec["source"], spec["filenames"], output_dir)

        return {"dir_name": output_dir}

    @flow
    def create_files():
        job1 = create_file("job1")
        job2 = create_file(
            "job2", copy=[{"source": job1["dir_name"], "filenames": "job1*"}]
        )
        return [job1, job2]

    scheduler.run(create_files())

    # Individual job folders/files should exist, and the job1 file should be
    # copied over to the job2 folder.
    assert Path(tmp_path / "job1/job1").exists()
    assert Path(tmp_path / "job2/job2").exists()
    assert Path(tmp_path / "job2/job1").exists()
