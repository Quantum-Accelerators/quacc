from __future__ import annotations

from pathlib import Path
from shutil import which

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose

from quacc.recipes.espresso.bands import bands_flow
from quacc.recipes.espresso.core import static_job
from quacc.utils.files import copy_decompress_files
from quacc.utils.kpts import convert_pmg_kpts

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("bands.x") is None, reason="QE not installed"
)

DATA_DIR = Path(__file__).parent / "data"


def test_bands_flow_conventional_cell(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, "Si.upf.gz", tmp_path)
    atoms = bulk("Si", cubic=True)
    calc_kwargs = {
        "input_data": {
            "control": {"pseudo_dir": tmp_path},
            "system": {"ecutwfc": 10, "ecutrho": 80},
            "electrons": {"conv_thr": 1.0e-6},
        },
        "pseudopotentials": {"Si": "Si.upf"},
    }
    scf_results = static_job(atoms, kpts=(2, 2, 2), **calc_kwargs)

    output = bands_flow(
        atoms,
        [{"source": scf_results["dir_name"], "filenames": "*"}],
        make_bandpath=True,
        line_density=1,
        job_params={"bands_pw_job": calc_kwargs},
    )

    bands_results = output["bands_pw"]
    assert len(bands_results["atoms"]) == 8
    assert_allclose(bands_results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-4)
    expected_kpts = convert_pmg_kpts({"line_density": 1}, atoms, force_gamma=True)[0]
    assert_allclose(bands_results["parameters"]["kpts"].cell, atoms.get_cell())
    assert_allclose(
        bands_results["parameters"]["kpts"].kpts, expected_kpts, atol=1.0e-8
    )
    assert Path(output["bands_pp"]["dir_name"], "bands.dat.gz").is_file()


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

    output = bands_flow(
        atoms,
        [{"source": tmp_path, "filenames": "*"}],
        line_density=1,
        job_params=job_params,
    )
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
        [{"source": tmp_path, "filenames": "*"}],
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
