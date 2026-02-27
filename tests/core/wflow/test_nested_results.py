from __future__ import annotations

import os
from pathlib import Path

import pytest
from ase.build import bulk

from quacc import change_settings, flow, get_settings, job, subflow
from quacc.wflow_tools.context import get_context_path, get_directory_context


def test_nested_bulk_slabs_results(tmp_path):
    """
    The final "auto-discovered" folder structure of a standard quacc recipe like `bulk_to_slabs_flow` that has a
    subflow/jobs looks like:

    <RESULTS_DIR>
    └── bulk_to_slabs_flow-<timestamp>
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
    with change_settings({"NESTED_RESULTS": True, "RESULTS_DIR": tmp_path}):
        # Make sure we import the `@flow` inside the context manager after
        # changing settings!
        from quacc.recipes.emt.slabs import bulk_to_slabs_flow

        atoms = bulk("Cu")
        bulk_to_slabs_flow(atoms, run_static=True)

        matches = list(
            tmp_path.glob(
                "bulk_to_slabs_flow*/bulk_to_slabs_subflow*/relax_job*/opt.json.gz"
            )
        )
        assert len(matches) == 4


def test_nested_results_job(tmp_path):
    with change_settings({"NESTED_RESULTS": True, "RESULTS_DIR": tmp_path}):

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
        assert any(tmp_path.glob("create_file_job*/foo"))


def test_nested_results_flow_subflows(tmp_path):
    with change_settings({"NESTED_RESULTS": True, "RESULTS_DIR": tmp_path}):

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
            tmp_path.glob("create_file_flow*/create_file_subflow*/create_file_job*/foo")
        )
        assert any(
            tmp_path.glob("create_file_flow*/create_file_subflow*/create_file_job*/bar")
        )


def test_nested_results_flow_no_subflows(tmp_path):
    with change_settings({"NESTED_RESULTS": True, "RESULTS_DIR": tmp_path}):

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
        assert any(tmp_path.glob("create_file_flow*/create_file_job*/foo"))
        assert any(tmp_path.glob("create_file_flow*/create_file_job*/bar"))


# Under no circumstances should `RESULTS_DIR` be nuked!
@pytest.mark.parametrize(
    ("nested_results", "scratch_dir", "create_unique_dir"),
    [
        (False, None, True),
        (True, None, True),
        (False, "use_tmp_path", True),
        (True, "use_tmp_path", True),
        (False, None, False),
        (True, None, False),
        (False, "use_tmp_path", False),
        (True, "use_tmp_path", False),
    ],
)
def test_results_dir_safe(tmp_path, nested_results, scratch_dir, create_unique_dir):
    # Replace the placeholder "use_tmp_path" with tmp_path
    scratch = tmp_path if scratch_dir == "use_tmp_path" else scratch_dir

    with change_settings(
        {
            "NESTED_RESULTS": nested_results,
            "RESULTS_DIR": tmp_path,
            "SCRATCH_DIR": scratch,
            "CREATE_UNIQUE_DIR": create_unique_dir,
        }
    ):
        # Make sure we import the `@flow` inside the context manager after
        # changing settings!
        from quacc.recipes.emt.slabs import bulk_to_slabs_flow

        atoms = bulk("Cu")
        bulk_to_slabs_flow(atoms, run_static=True)

        # The RESULTS_DIR should always exist
        assert tmp_path.exists()

        # No `tmp-` folders in RESULTS_DIR
        assert not any(x.name.startswith("tmp") for x in tmp_path.iterdir())


@pytest.mark.skipif(os.name == "nt", reason="Symlinks not supported on Windows")
@pytest.mark.parametrize("nested_results", [False, True])
def test_no_stale_symlinks_after_successful_flow(tmp_path, nested_results):
    """After a successful flow, all symlink-* in RESULTS_DIR should be removed."""
    scratch_dir = tmp_path / "scratch"
    results_dir = tmp_path / "results"
    scratch_dir.mkdir()
    results_dir.mkdir()

    with change_settings(
        {
            "NESTED_RESULTS": nested_results,
            "RESULTS_DIR": results_dir,
            "SCRATCH_DIR": scratch_dir,
        }
    ):
        from quacc.recipes.emt.slabs import bulk_to_slabs_flow

        atoms = bulk("Cu")
        bulk_to_slabs_flow(atoms, run_static=True)

        assert not any(x.name.startswith("symlink-") for x in results_dir.iterdir())


@pytest.mark.parametrize("nested_results", [False, True])
def test_no_symlinks_without_scratch_dir(tmp_path, nested_results):
    """When SCRATCH_DIR is not set, no symlinks should ever appear in RESULTS_DIR."""
    with change_settings(
        {"NESTED_RESULTS": nested_results, "RESULTS_DIR": tmp_path, "SCRATCH_DIR": None}
    ):
        from quacc.recipes.emt.slabs import bulk_to_slabs_flow

        atoms = bulk("Cu")
        bulk_to_slabs_flow(atoms, run_static=True)

        assert not any(x.is_symlink() for x in tmp_path.rglob("*"))


@pytest.mark.skipif(os.name == "nt", reason="Symlinks not supported on Windows")
@pytest.mark.parametrize("nested_results", [False, True])
def test_symlink_failed_after_job_failure(tmp_path, nested_results):
    """A failed job should leave a symlink-failed-* in RESULTS_DIR pointing to the failed dir."""
    from quacc import JobFailure

    scratch_dir = tmp_path / "scratch"
    results_dir = tmp_path / "results"
    scratch_dir.mkdir()
    results_dir.mkdir()

    with change_settings(
        {
            "NESTED_RESULTS": nested_results,
            "RESULTS_DIR": results_dir,
            "SCRATCH_DIR": scratch_dir,
        }
    ):
        from quacc.runners.prep import calc_setup, terminate

        _, _ = calc_setup(None)  # atoms=None to skip calc dir setup
        tmpdir = next(scratch_dir.glob("tmp-*"))  # grab the created tmpdir

        with pytest.raises(JobFailure):
            terminate(tmpdir, RuntimeError("failure"))

        failed_symlinks = list(results_dir.glob("symlink-failed-*"))
        assert len(failed_symlinks) == 1
        assert failed_symlinks[0].is_symlink()
        assert failed_symlinks[0].resolve().name.startswith("failed-")
        # Original tmp symlink should be gone
        assert not any(results_dir.glob("symlink-tmp-*"))
