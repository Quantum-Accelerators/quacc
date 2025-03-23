from __future__ import annotations

import pytest

parsl = pytest.importorskip("parsl")

from pathlib import Path

import numpy as np
from ase.build import bulk

from quacc import flow, job
from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412

TEST_RUNINFO = Path(__file__).parent.parent / "runinfo"
TEST_RENDER_DIR = Path(__file__).parent.parent / "_test_renders"


def test_render_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (4, 11, 7)
    atoms.rattle(stdev=0.5)
    atoms.set_cell(atoms.get_cell() * 1.2)

    result = relax_job(
        atoms,
        relax_cell=True,
        opt_params={
            "fmax": 0.01,  # Lower fmax and more steps for longer relaxation
            "max_steps": 500,
        },
        additional_fields={
            "render": {"output_dir": TEST_RENDER_DIR, "video_config": {"fps": 10}}
        },
    )

    assert result is not None
    # A trajectory mp4 should exist prefixed with Cu27 label
    assert any(
        f.name.startswith("Cu308") and f.suffix == ".mp4"
        for f in TEST_RENDER_DIR.glob("*.mp4")
    )


def test_render_relax_job_in_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    result = bulk_to_slabs_flow(
        atoms,
        job_params={
            "relax_job": {
                "opt_params": {"fmax": 0.1},
                "additional_fields": {"render": {"output_dir": TEST_RENDER_DIR}},
            }
        },
    )

    assert result is not None
    # A trajectory mp4 should exist prefixed with Cu96 label
    assert any(
        f.name.startswith("Cu96") and f.suffix == ".mp4"
        for f in TEST_RENDER_DIR.glob("*.mp4")
    )
