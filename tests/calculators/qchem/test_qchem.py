import os
from pathlib import Path

import pytest
from ase import units
from ase.io import read
from monty.io import zopen
from pymatgen.io.qchem.inputs import QCInput

from quacc.calculators.qchem import QChem

FILE_DIR = Path(__file__).resolve().parent
TEST_ATOMS = read(os.path.join(FILE_DIR, "test.xyz"))
OS_ATOMS = read(os.path.join(FILE_DIR, "OS_test.xyz"))


def test_qchem_write_input_basic(tmpdir):
    tmpdir.chdir()
    calc = QChem(TEST_ATOMS, cores=40)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 1
    calc.write_input(TEST_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "basic", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    assert not os.path.exists(os.path.join(FILE_DIR, "53.0"))

    with pytest.raises(NotImplementedError):
        QChem(TEST_ATOMS, cores=40, directory="notsupported")


def test_qchem_write_input_intermediate(tmpdir):
    tmpdir.chdir()
    params = {"dft_rung": 3, "basis_set": "def2-svpd", "pcm_dielectric": "3.0"}
    calc = QChem(TEST_ATOMS, cores=40, charge=-1, qchem_input_params=params)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == -1
    assert calc.parameters["spin_multiplicity"] == 2
    assert calc.parameters["dft_rung"] == 3
    assert calc.parameters["basis_set"] == "def2-svpd"
    assert calc.parameters["pcm_dielectric"] == "3.0"
    calc.write_input(TEST_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "intermediate", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()


def test_qchem_write_input_advanced(tmpdir):
    tmpdir.chdir()
    params = {
        "scf_algorithm": "gdm",
        "basis_set": "def2-svpd",
        "smd_solvent": "water",
        "overwrite_inputs": {"rem": {"method": "b97mv", "mem_total": "170000"}},
    }
    calc = QChem(
        TEST_ATOMS, cores=40, charge=-1, spin_multiplicity=2, qchem_input_params=params
    )
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == -1
    assert calc.parameters["spin_multiplicity"] == 2
    assert calc.parameters["scf_algorithm"] == "gdm"
    assert calc.parameters["basis_set"] == "def2-svpd"
    assert calc.parameters["smd_solvent"] == "water"
    assert calc.parameters["overwrite_rem_method"] == "b97mv"
    assert calc.parameters["overwrite_rem_mem_total"] == "170000"
    calc.write_input(TEST_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "advanced", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()


def test_qchem_write_input_open_shell_and_different_charges(tmpdir):
    tmpdir.chdir()
    calc = QChem(OS_ATOMS, cores=40)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 2
    calc.write_input(OS_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(os.path.join(FILE_DIR, "examples", "OSDC1.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(OS_ATOMS, cores=40, charge=0)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 2
    calc.write_input(OS_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(os.path.join(FILE_DIR, "examples", "OSDC1.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(OS_ATOMS, cores=40, charge=0, spin_multiplicity=4)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 0
    assert calc.parameters["spin_multiplicity"] == 4
    calc.write_input(OS_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(os.path.join(FILE_DIR, "examples", "OSDC2.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    with pytest.raises(ValueError):
        QChem(OS_ATOMS, spin_multiplicity=1)

    with pytest.raises(ValueError):
        QChem(OS_ATOMS, charge=0, spin_multiplicity=1)

    calc = QChem(OS_ATOMS, cores=40, charge=1)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 1
    assert calc.parameters["spin_multiplicity"] == 1
    calc.write_input(OS_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(os.path.join(FILE_DIR, "examples", "OSDC3.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(OS_ATOMS, cores=40, charge=1, spin_multiplicity=1)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 1
    assert calc.parameters["spin_multiplicity"] == 1
    calc.write_input(OS_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(os.path.join(FILE_DIR, "examples", "OSDC3.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()

    calc = QChem(OS_ATOMS, cores=40, charge=1, spin_multiplicity=3)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == 1
    assert calc.parameters["spin_multiplicity"] == 3
    calc.write_input(OS_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(os.path.join(FILE_DIR, "examples", "OSDC4.qin"))
    assert qcinp.as_dict() == ref_qcinp.as_dict()


def test_qchem_read_results_basic_and_write_53(tmpdir):
    tmpdir.chdir()
    calc = QChem(TEST_ATOMS, cores=40)
    os.chdir(os.path.join(FILE_DIR, "examples", "basic"))
    calc.read_results()
    assert calc.results["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-1.3826330655069403)
    assert calc.prev_orbital_coeffs is not None
    os.chdir(FILE_DIR)
    calc.write_input(TEST_ATOMS)
    assert os.path.exists(os.path.join(FILE_DIR, "53.0"))
    with zopen("53.0", mode="rb") as new_file:
        new_binary = new_file.read()
        with zopen(
            os.path.join(FILE_DIR, "examples", "basic", "53.0"), mode="rb"
        ) as old_file:
            old_binary = old_file.read()
            assert new_binary == old_binary
    qcinp = QCInput.from_file("mol.qin")
    assert qcinp.rem.get("scf_guess") == "read"


def test_qchem_read_results_intermediate(tmpdir):
    tmpdir.chdir()
    calc = QChem(TEST_ATOMS, cores=40)
    os.chdir(os.path.join(FILE_DIR, "examples", "intermediate"))
    calc.read_results()
    assert calc.results["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-0.6955571014353796)
    assert calc.prev_orbital_coeffs is not None


def test_qchem_read_results_advanced(tmpdir):
    tmpdir.chdir()
    calc = QChem(TEST_ATOMS, cores=40)
    os.chdir(os.path.join(FILE_DIR, "examples", "advanced"))
    calc.read_results()
    assert calc.results["energy"] == pytest.approx(-605.7310332390 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-0.4270884974249971)
    assert calc.prev_orbital_coeffs is not None
