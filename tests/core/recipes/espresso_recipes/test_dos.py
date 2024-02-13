from pathlib import Path
from shutil import which

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose

from quacc.recipes.espresso.dos import dos_flow, dos_job, projwfc_flow, projwfc_job
from quacc.utils.files import copy_decompress_files, copy_decompress_tree

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("dos.x") is None, reason="QE not installed"
)

DATA_DIR = Path(__file__).parent / "data"


def test_dos_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    copy_decompress_tree({DATA_DIR / "dos_test/": "pwscf.save/*.gz"}, tmp_path)
    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    output = dos_job(tmp_path)

    assert output["results"]["pwscf_dos"]["fermi"] == pytest.approx(7.199)


def test_projwfc_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    copy_decompress_tree({DATA_DIR / "dos_test/": "pwscf.save/*.gz"}, tmp_path)
    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    output = projwfc_job(tmp_path)
    print(output)
    assert output["name"] == "projwfc.x Projects-wavefunctions"
    assert output["parameters"]["input_data"]["projwfc"] == {}


def test_dos_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    atoms = bulk("Si")
    input_data = {
        "control": {"calculation": "scf", "pseudo_dir": tmp_path},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "static_job": {"input_data": input_data, "pseudopotentials": pseudopotentials},
        "non_scf_job": {"kspacing": 0.05},
    }

    output = dos_flow(atoms, job_params=job_params)
    assert_allclose(
        output["static_job"]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    assert (
        output["static_job"]["parameters"]["input_data"]["control"]["calculation"]
        == "scf"
    )
    assert (
        output["static_job"]["parameters"]["input_data"]["electrons"]["mixing_mode"]
        == "TF"
    )

    assert output["static_job"]["results"]["nbands"] == 8
    assert output["static_job"]["results"]["nspins"] == 1

    assert_allclose(
        output["non_scf_job"]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    assert (
        output["non_scf_job"]["parameters"]["input_data"]["control"]["calculation"]
        == "nscf"
    )
    assert (
        output["non_scf_job"]["parameters"]["input_data"]["electrons"]["mixing_mode"]
        == "TF"
    )
    assert output["non_scf_job"]["results"]["nbands"] == 8
    assert output["non_scf_job"]["results"]["nspins"] == 1

    assert output["dos_job"]["results"]["pwscf_dos"]["fermi"] == pytest.approx(6.772)


def test_projwfc_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    atoms = bulk("Si")
    input_data = {
        "control": {"calculation": "scf", "pseudo_dir": tmp_path},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "static_job": {"input_data": input_data, "pseudopotentials": pseudopotentials},
        "non_scf_job": {"kspacing": 0.05},
    }

    output = projwfc_flow(atoms, job_params=job_params)
    assert_allclose(
        output["static_job"]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    assert (
        output["static_job"]["parameters"]["input_data"]["control"]["calculation"]
        == "scf"
    )
    assert (
        output["static_job"]["parameters"]["input_data"]["electrons"]["mixing_mode"]
        == "TF"
    )

    assert output["static_job"]["results"]["nbands"] == 8
    assert output["static_job"]["results"]["nspins"] == 1

    assert_allclose(
        output["non_scf_job"]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    assert (
        output["non_scf_job"]["parameters"]["input_data"]["control"]["calculation"]
        == "nscf"
    )
    assert (
        output["non_scf_job"]["parameters"]["input_data"]["electrons"]["mixing_mode"]
        == "TF"
    )
    assert output["non_scf_job"]["results"]["nbands"] == 8
    assert output["non_scf_job"]["results"]["nspins"] == 1
