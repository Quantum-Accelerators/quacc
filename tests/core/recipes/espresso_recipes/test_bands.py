from __future__ import annotations

from pathlib import Path
from shutil import which

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose

from quacc.recipes.espresso.bands import bands_flow
from quacc.utils.files import copy_decompress_files

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("bands.x") is None, reason="QE not installed"
)

DATA_DIR = Path(__file__).parent / "data"


def test_bands_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR / "dos_test", Path("pwscf.save", "*.gz"), tmp_path)
    copy_decompress_files(DATA_DIR, "Si.upf.gz", tmp_path)
    atoms = bulk("Si")
    pseudopotentials = {"Si": "Si.upf"}
    job_params = {
        "bands_pw_job": {
            "input_data": {"control": {"pseudo_dir": tmp_path}},
            "pseudopotentials": pseudopotentials,
        },
        "fermi_surface_job": {"input_data": {"fermi": {}}},
    }

    output = bands_flow(atoms, tmp_path, line_density=1, job_params=job_params)
    assert (
        output["bands_pw"]["parameters"]["input_data"]["control"]["calculation"]
        == "bands"
    )

    assert Path(output["bands_pp"]["dir_name"], "bands.out.gz").exists()
    assert Path(output["bands_pp"]["dir_name"], "bands.dat.gz").exists()

    assert_allclose(
        output["bands_pw"]["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )

    assert output["bands_pw"]["results"]["nbands"] == 4

    assert output["bands_pp"]["name"] == "bands.x post-processing"


def test_bands_flow_with_fermi(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR / "dos_test", Path("pwscf.save", "*.gz"), tmp_path)
    copy_decompress_files(DATA_DIR, "Si.upf.gz", tmp_path)
    atoms = bulk("Si")
    pseudopotentials = {"Si": "Si.upf"}
    job_params = {
        "bands_pw_job": {
            "input_data": {"control": {"pseudo_dir": tmp_path}},
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.9,
        },
        "fermi_surface_job": {"input_data": {"fermi": {}}},
    }

    output = bands_flow(
        atoms,
        tmp_path,
        run_fermi_surface=True,
        make_bandpath=False,
        run_bands_pp=False,
        job_params=job_params,
    )
    assert (
        output["bands_pw"]["parameters"]["input_data"]["control"]["calculation"]
        == "bands"
    )

    assert_allclose(
        output["bands_pw"]["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )

    assert output["bands_pw"]["results"]["nbands"] == 4

    assert output["fermi_surface"]["name"] == "fs.x fermi_surface"
