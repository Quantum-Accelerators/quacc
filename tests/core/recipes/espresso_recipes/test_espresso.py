from pathlib import Path
from shutil import which

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose, assert_array_equal

from quacc import SETTINGS
from quacc.recipes.espresso.core import post_processing_job, relax_job, static_job
from quacc.recipes.espresso.phonons import phonon_job
from quacc.utils.files import copy_decompress_files

pytestmark = pytest.mark.skipif(which("pw.x") is None, reason="QE not installed")

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path}}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
    )

    assert_allclose(
        results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)
    assert_array_equal(
        results["atoms"].get_chemical_symbols(), atoms.get_chemical_symbols()
    )
    assert results["results"]["energy"] == pytest.approx(-293.71195934404255)

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["system"]["degauss"] == 0.001
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["system"]["ecutwfc"] == 30.0
    assert new_input_data["system"]["ecutrho"] == 240.0
    assert "kspacing" not in results["parameters"]
    assert results["parameters"].get("kpts") is None


def test_static_job_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

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

    assert_allclose(
        results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)
    assert_array_equal(
        results["atoms"].get_chemical_symbols(), atoms.get_chemical_symbols()
    )
    assert results["results"]["energy"] == pytest.approx(-293.71195934404255)

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"
    assert "kpts" not in results["parameters"]
    assert results["parameters"]["kspacing"] == 0.5

    pp_results = post_processing_job(prev_dir=results["dir_name"])
    assert Path(pp_results["dir_name"], "pseudo_charge_density.cube.gz").is_file()


def test_static_job_outdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {"pseudo_dir": tmp_path, "outdir": "test2/test3/test4"},
    }

    pseudopotentials = {"Si": "Si.upf"}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
    )

    assert_allclose(
        results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)
    assert_array_equal(
        results["atoms"].get_chemical_symbols(), atoms.get_chemical_symbols()
    )
    assert results["results"]["energy"] == pytest.approx(-293.71195934404255)

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"


def test_static_job_outdir_abs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

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
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
    )

    assert_allclose(
        results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)
    assert_array_equal(
        results["atoms"].get_chemical_symbols(), atoms.get_chemical_symbols()
    )
    assert results["results"]["energy"] == pytest.approx(-293.71195934404255)

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"


def test_static_job_dir_fail(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    with pytest.raises(ValueError):
        static_job(atoms, directory=Path("fake_path"))


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")
    atoms[0].position += 0.05

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path}}

    results = relax_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
    )

    with pytest.raises(AssertionError):
        assert_allclose(
            results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
        )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["control"]["calculation"] == "relax"


def test_relax_job_cell(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path}}

    results = relax_job(
        atoms,
        relax_cell=True,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=None,
    )

    with pytest.raises(AssertionError):
        assert_allclose(
            results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
        )
    with pytest.raises(AssertionError):
        assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["control"]["calculation"] == "vc-relax"


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

    sections = ["atoms", "eqpoints", "freqs", "kpoints", "mode_symmetries", "modes"]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]


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

    sections = ["atoms", "eqpoints", "freqs", "kpoints", "mode_symmetries", "modes"]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]
