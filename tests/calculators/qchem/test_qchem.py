import os
from pathlib import Path

import pytest
from ase import units
from ase.io import read
from pymatgen.io.qchem.inputs import QCInput

from quacc.calculators.qchem import QChem

FILE_DIR = Path(__file__).resolve().parent
TEST_ATOMS = read(os.path.join(FILE_DIR, "test.xyz"))


def test_qchem_write_input_basic():
    calc = QChem(TEST_ATOMS, 40)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] is None
    assert calc.parameters["spin_multiplicity"] is None
    calc.write_input(TEST_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "basic", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    os.remove("mol.qin")


def test_qchem_write_input_intermediate():
    params = {"dft_rung": 3, "basis_set": "def2-svpd", "pcm_dielectric": "3.0"}
    calc = QChem(TEST_ATOMS, cores=40, charge=-1, qchem_input_params=params)
    assert calc.parameters["cores"] == 40
    assert calc.parameters["charge"] == -1
    assert calc.parameters["spin_multiplicity"] is None
    assert calc.parameters["dft_rung"] == 3
    assert calc.parameters["basis_set"] == "def2-svpd"
    assert calc.parameters["pcm_dielectric"] == "3.0"
    calc.write_input(TEST_ATOMS)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "intermediate", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    os.remove("mol.qin")


def test_qchem_write_input_advanced():
    params = {
        "scf_algorithm": "gdm",
        "qchem_version": 6,
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
    assert calc.parameters["qchem_version"] == 6
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
    os.remove("mol.qin")


def test_qchem_read_results_basic():
    calc = QChem(TEST_ATOMS, 40)
    os.chdir(os.path.join(FILE_DIR, "examples", "basic"))
    calc.read_results()
    assert calc.results["energy"] == pytest.approx(-606.1616819641 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-1.3826330655069403)


def test_qchem_read_results_intermediate():
    calc = QChem(TEST_ATOMS, 40)
    os.chdir(os.path.join(FILE_DIR, "examples", "intermediate"))
    calc.read_results()
    assert calc.results["energy"] == pytest.approx(-605.6859554025 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-0.6955571014353796)


def test_qchem_read_results_advanced():
    calc = QChem(TEST_ATOMS, 40)
    os.chdir(os.path.join(FILE_DIR, "examples", "advanced"))
    calc.read_results()
    assert calc.results["energy"] == pytest.approx(-605.7310332390 * units.Hartree)
    assert calc.results["forces"][0][0] == pytest.approx(-0.4270884974249971)
