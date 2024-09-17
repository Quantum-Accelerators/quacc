from __future__ import annotations

import os
from importlib.util import find_spec
from pathlib import Path

import pytest
from ase import units
from ase.io import read
from pymatgen.io.qchem.inputs import QCInput

from quacc.calculators.qchem import QChem

has_obabel = bool(find_spec("openbabel"))

FILE_DIR = Path(__file__).parent


@pytest.fixture
def test_atoms():
    return read(FILE_DIR / "test.xyz")


@pytest.fixture
def os_atoms():
    return read(FILE_DIR / "OS_test.xyz")


def test_qchem_write_input_basic(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(
        test_atoms,
        qchem_dict_set_params={
            "basis_set": "def2-tzvpd",
            "job_type": "force",
            "scf_algorithm": "diis",
        },
    )
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 1
    calc.write_input(test_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "basic" / "mol.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    assert not Path(FILE_DIR / "53.0").exists()

    with pytest.raises(
        NotImplementedError,
        match="Do not specify `molecule` in `qchem_dict_set_params`",
    ):
        QChem(
            test_atoms, job_type="freq", qchem_dict_set_params={"molecule": "test"}
        ).write_input(test_atoms)


def test_qchem_write_input_intermediate(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(
        test_atoms,
        qchem_dict_set_params={
            "job_type": "force",
            "basis_set": "def2-svpd",
            "dft_rung": 3,
            "pcm_dielectric": "3.0",
            "scf_algorithm": "diis",
        },
        charge=-1,
        spin_multiplicity=2,
    )
    assert calc.parameters["charge"] == -1
    assert calc.parameters["spin_multiplicity"] == 2
    calc.write_input(test_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        str(FILE_DIR / "examples" / "intermediate" / "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()


def test_qchem_write_input_advanced(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(
        test_atoms,
        qchem_dict_set_params={
            "job_type": "force",
            "basis_set": "def2-svpd",
            "scf_algorithm": "gdm",
            "max_scf_cycles": 100,
            "smd_solvent": "water",
            "overwrite_inputs": {"rem": {"method": "b97mv", "mem_total": "170000"}},
        },
        charge=-1,
        spin_multiplicity=2,
    )
    assert calc.parameters["charge"] == -1
    assert calc.parameters["spin_multiplicity"] == 2
    assert "method" not in calc.parameters
    calc.write_input(test_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "advanced" / "mol.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()


def test_qchem_write_input_open_shell_and_different_charges(
    tmp_path, monkeypatch, os_atoms
):
    monkeypatch.chdir(tmp_path)
    calc = QChem(
        os_atoms,
        spin_multiplicity=2,
        qchem_dict_set_params={
            "basis_set": "def2-tzvpd",
            "job_type": "force",
            "scf_algorithm": "diis",
        },
    )
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 2
    calc.write_input(os_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC1.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(
        os_atoms,
        charge=0,
        spin_multiplicity=4,
        qchem_dict_set_params={
            "basis_set": "def2-tzvpd",
            "job_type": "force",
            "scf_algorithm": "diis",
        },
    )
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 4
    calc.write_input(os_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC2.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(
        os_atoms,
        charge=1,
        qchem_dict_set_params={
            "basis_set": "def2-tzvpd",
            "job_type": "force",
            "scf_algorithm": "diis",
        },
    )
    assert calc.parameters["charge"] == 1
    assert calc.parameters["spin_multiplicity"] == 1
    calc.write_input(os_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC3.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(
        os_atoms,
        charge=1,
        spin_multiplicity=3,
        qchem_dict_set_params={
            "basis_set": "def2-tzvpd",
            "job_type": "force",
            "scf_algorithm": "diis",
        },
    )
    assert calc.parameters["charge"] == 1
    assert calc.parameters["spin_multiplicity"] == 3
    calc.write_input(os_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "OSDC4.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()


def test_qchem_write_input_freq(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(
        test_atoms,
        qchem_dict_set_params={
            "scf_algorithm": "diis",
            "job_type": "freq",
            "basis_set": "def2-svpd",
            "dft_rung": 3,
            "pcm_dielectric": 3.0,
        },
        charge=-1,
        spin_multiplicity=2,
    )
    assert calc.parameters["charge"] == -1
    assert calc.parameters["spin_multiplicity"] == 2
    calc.write_input(test_atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(str(FILE_DIR / "examples" / "freq" / "mol.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()


@pytest.mark.skipif(has_obabel is False, reason="openbabel needed")
def test_qchem_read_results_basic_and_write_53(tmp_path, monkeypatch, test_atoms):
    calc = QChem(
        test_atoms,
        rem={"basis": "def2-tzvpd", "method": "wb97x-v", "job_type": "force"},
    )
    monkeypatch.chdir(FILE_DIR / "examples" / "basic")
    calc.read_results()
    monkeypatch.chdir(tmp_path)

    assert calc.results["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-1.3826330655069403)
    assert calc.prev_orbital_coeffs is not None

    calc.write_input(test_atoms)
    assert Path(tmp_path, "53.0").exists()
    with open("53.0", mode="rb") as new_file:
        new_binary = new_file.read()
        with open(
            os.path.join(FILE_DIR, "examples", "basic", "53.0"), mode="rb"
        ) as old_file:
            old_binary = old_file.read()
            assert new_binary == old_binary
    qcinp = QCInput.from_file("mol.qin")
    assert qcinp.rem.get("scf_guess") == "read"


@pytest.mark.skipif(has_obabel is False, reason="openbabel needed")
def test_qchem_read_results_intermediate(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(test_atoms)
    monkeypatch.chdir(FILE_DIR / "examples" / "intermediate")
    calc.read_results()

    assert calc.results["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-0.6955571014353796)
    assert calc.prev_orbital_coeffs is not None


@pytest.mark.skipif(has_obabel is False, reason="openbabel needed")
def test_qchem_read_results_advanced(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(test_atoms)
    monkeypatch.chdir(FILE_DIR / "examples" / "advanced")
    calc.read_results()

    assert calc.results["energy"] == pytest.approx(-605.7310332390 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-0.4270884974249971)
    assert calc.prev_orbital_coeffs is not None
    assert calc.results.get("hessian") is None


@pytest.mark.skipif(has_obabel is False, reason="openbabel needed")
def test_qchem_read_results_freq(tmp_path, monkeypatch, test_atoms):
    monkeypatch.chdir(tmp_path)
    calc = QChem(test_atoms, job_type="freq")
    monkeypatch.chdir(FILE_DIR / "examples" / "freq")
    calc.read_results()

    assert calc.results["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert calc.results.get("forces") is None
    assert calc.prev_orbital_coeffs is not None
    assert len(calc.results["hessian"]) == 42
    assert len(calc.results["hessian"][0]) == 42
    assert calc.results["taskdoc"]["output"]["frequencies"][0] == -340.2
    assert len(calc.results["taskdoc"]["output"]["frequencies"]) == 36
    assert len(calc.results["taskdoc"]["output"]["frequency_modes"]) == 36
    assert len(calc.results["taskdoc"]["output"]["frequency_modes"][0]) == 14
    assert len(calc.results["taskdoc"]["output"]["frequency_modes"][0][0]) == 3
