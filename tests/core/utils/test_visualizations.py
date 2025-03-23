from __future__ import annotations

import pytest

parsl = pytest.importorskip("parsl")

from pathlib import Path

from ase.build import bulk

from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412

TEST_RUNINFO = Path(__file__).parent.parent / "runinfo"
TEST_RENDER_DIR = Path(__file__).parent.parent / "_test_renders"


def test_render_option_variations(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    # Test with None render config
    result_no_render = relax_job(
        atoms,
        relax_cell=True,
        opt_params={"fmax": 0.8},
        additional_fields={"render": None},
    )
    assert result_no_render is not None

    # Test with empty render config (no output_dir)
    result_no_output = relax_job(
        atoms,
        relax_cell=True,
        opt_params={"fmax": 0.8},
        additional_fields={"render": {}},
    )
    assert result_no_output is not None

    # Test with no render config at all
    result_no_config = relax_job(
        atoms, relax_cell=True, opt_params={"fmax": 0.8}, additional_fields={}
    )
    assert result_no_config is not None


def test_render_within_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 7, 5)
    atoms.rattle(stdev=0.5)
    atoms.set_cell(atoms.get_cell() * 1.2)

    result = relax_job(
        atoms,
        relax_cell=True,
        opt_params={"fmax": 0.8},
        additional_fields={
            "render": {"output_dir": TEST_RENDER_DIR, "video_config": {"fps": 10}}
        },
    )

    assert result is not None
    # A trajectory mp4 should exist prefixed with Cu27 label
    assert any(
        f.name.startswith("Cu70") and f.suffix == ".mp4"
        for f in TEST_RENDER_DIR.glob("*.mp4")
    )
    # A final image should exist prefixed with Cu96 label
    assert any(
        f.name.startswith("Cu70") and f.suffix == ".png"
        for f in TEST_RENDER_DIR.glob("*.png")
    )


def test_render_within_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    result = bulk_to_slabs_flow(
        atoms,
        job_params={
            "relax_job": {
                "opt_params": {"fmax": 0.1},
                "additional_fields": {
                    "render": {
                        "output_dir": TEST_RENDER_DIR,
                        "video_config": {
                            "width": 801,
                            "height": 603,
                        },  # ensuring odd dimensions don't break ffmpeg
                    }
                },
            }
        },
    )

    assert result is not None
    # A trajectory mp4 should exist prefixed with Cu96 label
    assert any(
        f.name.startswith("Cu96") and f.suffix == ".mp4"
        for f in TEST_RENDER_DIR.glob("*.mp4")
    )
