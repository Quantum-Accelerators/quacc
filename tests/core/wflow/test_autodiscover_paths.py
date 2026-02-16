from __future__ import annotations

from pathlib import Path

from ase.build import bulk

from quacc import change_settings, flow, get_settings, job, subflow
from quacc.recipes.emt.slabs import bulk_to_slabs_flow
from quacc.wflow_tools.context import get_context_path, get_directory_context


def test_autodiscover_bulk_slabs_paths(tmp_path):
    """
    The final "auto-discovered" folder structure of a standard quacc recipe like `bulk_to_slabs_flow` that has a
    subflow/jobs looks like:

    <RESULTS_DIR>
    └── quacc-<timestamp>
        └── bulk_to_slabs_flow
            └── bulk_to_slabs_subflow-<timestamp>
                ├── relax_job-<timestamp>
                │         ├── opt.json.gz
                │         ├── opt.log.gz
                │         └── opt.traj.gz
                ├── relax_job-<timestamp>
                │         ├── opt.json.gz
                │         ├── opt.log.gz
                │         └── opt.traj.gz
                ...
                ├── static_job-<timestamp>
                ...
                └── static_job-<timestamp>
    """
    with change_settings({"RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu")
        bulk_to_slabs_flow(atoms, run_static=True)

        matches = list(
            tmp_path.glob(
                "quacc-*/bulk_to_slabs_flow/bulk_to_slabs_subflow*/relax_job*/opt.json.gz"
            )
        )
        assert len(matches) == 4


def test_autodiscover_paths_job(tmp_path):
    with change_settings({"RESULTS_DIR": tmp_path}):

        @job
        def create_file_job(name):
            settings = get_settings()
            job_results_dir = (
                settings.RESULTS_DIR
                / Path(get_directory_context())
                / get_context_path()
            )
            job_results_dir.mkdir(parents=True, exist_ok=True)
            job_results_file = job_results_dir / f"{name}"
            job_results_file.touch()

        create_file_job("foo")
        assert any(tmp_path.glob("quacc-*/create_file_job/foo"))


def test_autodiscover_paths_flow_subflows(tmp_path):
    with change_settings({"RESULTS_DIR": tmp_path}):

        @job
        def create_file_job(name):
            settings = get_settings()
            job_results_dir = (
                settings.RESULTS_DIR
                / Path(get_directory_context())
                / get_context_path()
            )
            job_results_dir.mkdir(parents=True, exist_ok=True)
            job_results_file = job_results_dir / f"{name}"
            job_results_file.touch()

        @subflow
        def create_file_subflow():
            return [create_file_job("foo"), create_file_job("bar")]

        @flow
        def create_file_flow():
            # Return 2 identical subflows
            return [create_file_subflow(), create_file_subflow()]

        create_file_flow()
        assert any(
            tmp_path.glob(
                "quacc-*/create_file_flow/create_file_subflow*/create_file_job*/foo"
            )
        )
        assert any(
            tmp_path.glob(
                "quacc-*/create_file_flow/create_file_subflow*/create_file_job*/bar"
            )
        )


def test_autodiscover_paths_flow_no_subflows(tmp_path):
    with change_settings({"RESULTS_DIR": tmp_path}):

        @job
        def create_file_job(name):
            settings = get_settings()
            job_results_dir = (
                settings.RESULTS_DIR
                / Path(get_directory_context())
                / get_context_path()
            )
            job_results_dir.mkdir(parents=True, exist_ok=True)
            job_results_file = job_results_dir / f"{name}"
            job_results_file.touch()

        @flow
        def create_file_flow():
            return [create_file_job("foo"), create_file_job("bar")]

        create_file_flow()
        assert any(tmp_path.glob("quacc-*/create_file_flow/create_file_job*/foo"))
        assert any(tmp_path.glob("quacc-*/create_file_flow/create_file_job*/bar"))


def test_autodiscover_paths_flow_replace(tmp_path):
    with change_settings({"RESULTS_DIR": tmp_path}):

        @job
        def create_file_job_replaced(name):
            settings = get_settings()
            job_results_dir = (
                settings.RESULTS_DIR
                / Path(get_directory_context())
                / get_context_path()
            )
            job_results_dir.mkdir(parents=True, exist_ok=True)
            job_results_file = job_results_dir / f"{name}"
            job_results_file.touch()

        @job
        def create_file_job(name):
            # When a job replaces another, a level is added to the folder hierarchy
            return create_file_job_replaced(name)

        @flow
        def create_file_flow():
            return create_file_job("foo")

        create_file_flow()
        assert any(
            tmp_path.glob(
                "quacc-*/create_file_flow/create_file_job*/create_file_job_replaced*/foo"
            )
        )
