from __future__ import annotations

import pytest

jf = pytest.importorskip("jobflow")

from pathlib import Path

from ase.build import bulk

from quacc import change_settings
from quacc.wflow_tools.context import get_context_path, get_directory_context


def test_autodiscover_bulk_slabs_paths(tmp_path):
    atoms = bulk("Cu")

    with change_settings({"AUTODISCOVER_DIR": True, "RESULTS_DIR": tmp_path}):
        # Make sure we import the `@flow` inside the context manager after
        # changing settings!
        from quacc.recipes.emt.slabs import bulk_to_slabs_flow

        flow = bulk_to_slabs_flow(atoms, run_static=False)
        jf.run_locally(flow, create_folders=False, ensure_success=True)

        matches = list(
            tmp_path.glob(
                "quacc-*/bulk_to_slabs_flow/bulk_to_slabs_subflow/relax_job-*/opt.log.gz"
            )
        )
        assert len(matches) == 4


def test_autodiscover_paths(tmp_path):
    with change_settings({"AUTODISCOVER_DIR": True, "RESULTS_DIR": tmp_path}):
        from quacc import flow, job

        @job
        def create_file_job(name):
            job_results_dir = Path(get_directory_context()) / get_context_path()
            job_results_dir.mkdir(parents=True, exist_ok=True)
            job_results_file = job_results_dir / f"{name}"
            job_results_file.touch()

        @flow
        def create_file_flow():
            return [create_file_job("foo"), create_file_job("bar")]

        flow = create_file_flow()
        jf.run_locally(flow, create_folders=False, ensure_success=True)
        assert any(tmp_path.glob("quacc-*/create_file_flow/create_file_job/foo"))
        assert any(tmp_path.glob("quacc-*/create_file_flow/create_file_job/bar"))


def test_autodiscover_paths_flow_replace(tmp_path):
    with change_settings({"AUTODISCOVER_DIR": True, "RESULTS_DIR": tmp_path}):
        from quacc import flow, job

        @job
        def create_file_job_replaced(name):
            job_results_dir = Path(get_directory_context()) / get_context_path()
            job_results_dir.mkdir(parents=True, exist_ok=True)
            job_results_file = job_results_dir / f"{name}"
            job_results_file.touch()

        @job
        def create_file_job(name):
            return create_file_job_replaced(name)

        @flow
        def create_file_flow():
            return create_file_job("foo")

        flow = create_file_flow()
        jf.run_locally(flow, create_folders=False, ensure_success=True)
        # Replaced job(s) show up a subfolder(s) within the original job folder
        assert any(
            tmp_path.glob(
                "quacc-*/create_file_flow/create_file_job/create_file_job_replaced/foo"
            )
        )
