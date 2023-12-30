from pathlib import Path
from shutil import which

import numpy as np
import pytest
from ase.build import bulk
from ase.io.espresso import construct_namelist

from quacc import SETTINGS
from quacc.recipes.espresso.core import static_job
from quacc.recipes.espresso.phonons import phonon_job, grid_phonon_flow
from quacc.utils.files import copy_decompress_files

pytestmark = pytest.mark.skipif(which("pw.x") is None, reason="QE not installed")

DEFAULT_SETTINGS = SETTINGS.model_copy()


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path}}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    assert np.allclose(results["atoms"].positions, atoms.positions, atol=1.0e-4)

    assert np.allclose(results["atoms"].cell, atoms.cell, atol=1.0e-3)
    assert (results["atoms"].symbols == atoms.symbols).all()

    new_input_data = results["parameters"]["input_data"]

    assert new_input_data["system"]["degauss"] == 0.001
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["system"]["ecutwfc"] == 30.0
    assert new_input_data["system"]["ecutrho"] == 240.0


def test_static_job_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {"pseudo_dir": tmp_path},
    }

    pseudopotentials = {"Si": "Si.upf"}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    assert np.allclose(results["atoms"].positions, atoms.positions, atol=1.0e-4)

    assert np.allclose(results["atoms"].cell, atoms.cell, atol=1.0e-3)
    assert (results["atoms"].symbols == atoms.symbols).all()

    new_input_data = results["parameters"]["input_data"]

    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"


def test_static_job_outdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {"pseudo_dir": tmp_path, "outdir": "test2/test3/test4"},
    }

    pseudopotentials = {"Si": "Si.upf"}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    input_data = dict(construct_namelist(input_data))

    assert np.allclose(results["atoms"].positions, atoms.positions, atol=1.0e-4)

    assert np.allclose(results["atoms"].cell, atoms.cell, atol=1.0e-3)
    assert (results["atoms"].symbols == atoms.symbols).all()

    new_input_data = results["parameters"]["input_data"]

    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"


def test_static_job_outdir_abs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {
            "pseudo_dir": tmp_path,
            "outdir": Path(tmp_path, "test2").resolve(),
        },
    }
    pseudopotentials = {"Si": "Si.upf"}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    input_data = dict(construct_namelist(input_data))

    assert np.allclose(results["atoms"].positions, atoms.positions, atol=1.0e-4)

    assert np.allclose(results["atoms"].cell, atoms.cell, atol=1.0e-3)
    assert (results["atoms"].symbols == atoms.symbols).all()

    new_input_data = results["parameters"]["input_data"]

    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"


def test_static_job_dir_fail(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    with pytest.raises(ValueError):
        static_job(atoms, directory=Path("fake_path"))


def test_phonon_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-8}}

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.25
    )

    ph_results = phonon_job(pw_results["dir_name"], input_data=ph_loose)

    assert (0, 0, 0) in ph_results["results"]
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].positions,
        atoms.positions,
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].cell, atoms.cell, atol=1.0e-3
    )
    assert (ph_results["results"][(0, 0, 0)]["atoms"].symbols == atoms.symbols).all()

    sections = ["atoms", "eqpoints", "freqs", "kpoints", "mode_symmetries", "modes"]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]


def test_phonon_job_list_to_do(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    ph_loose = {
        "inputph": {"tr2_ph": 1e-8, "qplot": True, "nat_todo": 1, "ldisp": True}
    }

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.25
    )

    qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

    nat_todo = [1]

    ph_results = phonon_job(
        pw_results["dir_name"], input_data=ph_loose, qpts=qpts, nat_todo=nat_todo
    )

    assert (0, 0, 0) in ph_results["results"]
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].positions,
        atoms.positions,
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].cell, atoms.cell, atol=1.0e-3
    )
    assert (ph_results["results"][(0, 0, 0)]["atoms"].symbols == atoms.symbols).all()

    sections = ["atoms", "eqpoints", "freqs", "kpoints", "mode_symmetries", "modes"]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]

def test_phonon_grid(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-8}}

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {'pw_job': {'input_data' : input_data, 'pseudopotentials' : pseudopotentials}, 'ph_job': {'input_data': ph_loose}}

    grid_results = grid_phonon_flow(atoms, job_params = job_params)
