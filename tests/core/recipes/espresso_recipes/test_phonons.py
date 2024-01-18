from pathlib import Path
from shutil import which

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose, assert_array_equal

from quacc import SETTINGS
from quacc.recipes.espresso.core import static_job
from quacc.recipes.espresso.phonons import grid_phonon_flow, phonon_job
from quacc.utils.files import copy_decompress_files

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("ph.x") is None, reason="QE not installed"
)

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_phonon_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    copy_decompress_files([DATA_DIR / "Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-8}}

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    ph_results = phonon_job(pw_results["dir_name"], input_data=ph_loose)

    assert (0, 0, 0) in ph_results["results"]
    assert_allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert_allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].get_cell(),
        atoms.get_cell(),
        atol=1.0e-3,
    )
    assert_array_equal(
        ph_results["results"][(0, 0, 0)]["atoms"].get_chemical_symbols(),
        atoms.get_chemical_symbols(),
    )

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_job_list_to_do(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    copy_decompress_files([DATA_DIR / "Li.upf.gz"], tmp_path)

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
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

    nat_todo = [1]

    ph_results = phonon_job(
        pw_results["dir_name"],
        input_data=ph_loose,
        qpts=qpts,
        nat_todo_indices=nat_todo,
    )

    assert (0, 0, 0) in ph_results["results"]
    assert_allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert_allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].get_cell(),
        atoms.get_cell(),
        atol=1.0e-3,
    )
    assert_array_equal(
        ph_results["results"][(0, 0, 0)]["atoms"].get_chemical_symbols(),
        atoms.get_chemical_symbols(),
    )

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_grid(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-6, "qplot": True, "ldisp": True}}

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose, "qpts": [(0, 0, 0, 1), (0.5, 0.0, 0.0, 1)]},
    }

    grid_results = grid_phonon_flow(atoms, job_params=job_params)

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][(0, 0, 0)]


def test_phonon_grid_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Li.upf.gz"], tmp_path)

    atoms = bulk("Li", "bcc", orthorhombic=True)

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-6}}

    pseudopotentials = {"Li": "Li.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose},
    }

    grid_results = grid_phonon_flow(atoms, job_params=job_params, nblocks=3)

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][(0, 0, 0)]
