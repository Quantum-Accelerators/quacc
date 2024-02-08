from shutil import which

import pytest

from quacc import SETTINGS

pytestmark = pytest.mark.skipif(
    which(str(SETTINGS.ESPRESSO_BINARIES["pw"])) is None, reason="QE not installed"
)

from pathlib import Path
from shutil import which
from subprocess import CalledProcessError

import pytest
from ase.build import bulk
from ase.optimize import BFGS
from numpy.testing import assert_allclose, assert_array_equal

from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.calculators.espresso.utils import pw_copy_files
from quacc.recipes.espresso.core import (
    ase_relax_job,
    non_scf_job,
    post_processing_job,
    relax_job,
    static_job,
)
from quacc.utils.files import copy_decompress_files

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

    with pytest.raises(ValueError):
        static_job(
            atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
        )


def test_static_job_dir_fail(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    with pytest.raises(ValueError):
        static_job(atoms, directory=Path("fake_path"))


def test_static_job_test_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}

    EspressoTemplate._test_run({"input_data": {"control": {"prefix": "test"}}}, Path())

    assert Path("test.EXIT").exists()

    with pytest.raises(CalledProcessError):
        static_job(
            atoms,
            pseudopotentials=pseudopotentials,
            input_data={"pseudo_dir": tmp_path},
            test_run=True,
        )


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


def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")
    atoms[0].position += 0.05

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path}}

    results = ase_relax_job(
        atoms,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=None,
        opt_params={"max_steps": 10, "fmax": 1.0e-1, "optimizer": BFGS},
    )

    with pytest.raises(AssertionError):
        assert_allclose(
            results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
        )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)

    assert len(results["trajectory"]) == 3
    assert_allclose(
        results["trajectory_results"][-1]["free_energy"], -293.71198070497894
    )
    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["control"]["calculation"] == "scf"


def test_ase_relax_cell_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")
    atoms[0].position += 0.05

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {
        "control": {"pseudo_dir": tmp_path},
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.005,
    }

    with pytest.raises(RuntimeError):
        ase_relax_job(
            atoms,
            relax_cell=True,
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kpts=None,
            opt_params={"max_steps": 2, "fmax": 1.0e-1, "optimizer": BFGS},
        )


def test_relax_job_cell(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {
        "control": {"pseudo_dir": tmp_path, "etot_conv_thr": 1.0},
        "system": {"occupations": "smearing", "degauss": 0.001},
        "cell": {"press_conv_thr": 1.5e2},
    }

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


def test_non_scf_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path}}
    static_result = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
    )
    file_to_copy = pw_copy_files(
        input_data, static_result["dir_name"], include_wfc=False
    )
    results = non_scf_job(
        atoms, file_to_copy, input_data=input_data, pseudopotentials=pseudopotentials
    )

    assert_allclose(
        results["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )
    assert_allclose(results["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3)
    assert_array_equal(
        results["atoms"].get_chemical_symbols(), atoms.get_chemical_symbols()
    )
    assert results["parameters"]["input_data"]["control"]["calculation"] == "nscf"
    assert results["results"]["nspins"] == 1
    assert results["results"]["nbands"] == 4

    new_input_data = results["parameters"]["input_data"]
    assert new_input_data["system"]["ecutwfc"] == 30.0
    assert new_input_data["system"]["ecutrho"] == 240.0
    assert "kspacing" not in results["parameters"]
    assert results["parameters"].get("kpts") is None


def test_pw_copy(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    pseudopotentials = {"Si": "Si.upf"}
    input_data = {"control": {"pseudo_dir": tmp_path, "max_seconds": 5}}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=None
    )

    new_input_data = results["parameters"]["input_data"]
    new_input_data["restart_mode"] = "restart"
    new_input_data["max_seconds"] = 10**7

    files_to_copy = pw_copy_files(new_input_data, results["dir_name"], include_wfc=True)

    results = static_job(
        atoms,
        input_data=new_input_data,
        pseudopotentials=pseudopotentials,
        kpts=None,
        copy_files=files_to_copy,
    )
    assert new_input_data["system"]["ecutwfc"] == 30.0
    assert new_input_data["system"]["ecutrho"] == 240.0
