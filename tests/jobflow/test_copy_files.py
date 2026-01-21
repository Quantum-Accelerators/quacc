from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")

from pathlib import Path

from quacc import job
from quacc.wflow_tools.job_argument import Copy
from quacc.utils.files import copy_decompress_files


def test_copy_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def create_file(name: str, copy: Copy | None = None):
        output_dir = tmp_path / name
        output_dir.mkdir(parents=True, exist_ok=True)
        Path(output_dir / name).touch()

        if copy is not None:
            uuids_to_src_dirs = copy.get("uuids_to_keys", {})
            uuids_to_files = copy.get("uuids_to_vals", {})
            for src_dir_uuid, src_dir in uuids_to_src_dirs.items():
                copy_decompress_files(src_dir, uuids_to_files.get(src_dir_uuid, []), output_dir)

        return {"dir_name": output_dir}

    job1 = create_file("job1")
    job2 = create_file(
        "job2", copy=Copy({job1.output["dir_name"]: "job1*"})
    )
    flow = jf.Flow([job1, job2])
    assert "DiGraph with 2 nodes and 1 edges" == str(flow.graph)
    jf.run_locally(flow, ensure_success=True, create_folders=True)

    # Individual job folders/files should exist, and the job1 file should be
    # copied over to the job2 folder.
    assert Path(tmp_path / "job1/job1").exists()
    assert Path(tmp_path / "job2/job2").exists()
    assert Path(tmp_path / "job2/job1").exists()
