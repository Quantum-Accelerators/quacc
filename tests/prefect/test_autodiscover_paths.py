from __future__ import annotations

import os
from pathlib import Path

from ase.build import bulk

from quacc import flow, job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow
from quacc.wflow_tools.context import get_context_path, get_directory_context


def test_autodiscover_bulk_slabs_paths():
    atoms = bulk("Cu")
    bulk_to_slabs_flow(atoms)

    results_dir = os.environ["QUACC_RESULTS_DIR"]

    matches = list(
        Path(results_dir).glob(
            "quacc-*/bulk_to_slabs_flow/bulk_to_slabs_subflow/relax_job-*/opt.json.gz"
        )
    )
    assert len(matches) == 4

    matches = list(
        Path(results_dir).glob(
            "quacc-*/bulk_to_slabs_flow/bulk_to_slabs_subflow/static_job-*"
        )
    )
    assert len(matches) == 4


def test_autodiscover_paths():
    @job
    def create_file_job(name):
        job_results_dir = Path(get_directory_context()) / get_context_path()
        job_results_dir.mkdir(parents=True, exist_ok=True)
        job_results_file = job_results_dir / f"{name}"
        job_results_file.touch()

    @flow
    def create_file_flow():
        return [create_file_job("foo"), create_file_job("bar")]

    create_file_flow()

    results_dir = os.environ["QUACC_RESULTS_DIR"]
    assert any(Path(results_dir).glob("quacc-*/create_file_flow/create_file_job/foo"))
    assert any(Path(results_dir).glob("quacc-*/create_file_flow/create_file_job/bar"))
