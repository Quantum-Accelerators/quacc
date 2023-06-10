import os
from copy import deepcopy
from pathlib import Path

import ase.units
import numpy as np
import pytest
from ase.io import read
from pymatgen.core import Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

from quacc.calculators.qchem import QChem
from quacc.util.atoms import prep_next_run

FILE_DIR = Path(__file__).resolve().parent


def test_qchem_write_input_basic():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    calc = QChem(atoms)
    calc.write_input(atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "basic", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    os.remove("mol.qin")


def test_qchem_write_input_intermediate():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    params = {"dft_rung": 3, "basis_set": "def2-svpd", "pcm_dielectric": "3.0"}
    calc = QChem(atoms, charge=-1, qchem_input_params=params)
    calc.write_input(atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "intermediate", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    os.remove("mol.qin")


def test_qchem_write_input_advanced():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    params = {
        "scf_algorithm": "gdm",
        "qchem_version": 6,
        "dft_rung": 3,
        "basis_set": "def2-svpd",
        "smd_solvent": "water",
        "overwrite_inputs": {"rem": {"mem_total": "170000"}},
    }
    calc = QChem(atoms, charge=-1, spin_multiplicity=2, qchem_input_params=params)
    calc.write_input(atoms)
    qcinp = QCInput.from_file("mol.qin")
    ref_qcinp = QCInput.from_file(
        os.path.join(FILE_DIR, "examples", "advanced", "mol.qin")
    )
    assert qcinp.as_dict() == ref_qcinp.as_dict()
    os.remove("mol.qin")


def test_qchem_read_results_basic():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    calc = QChem(atoms)
    os.chdir(os.path.join(FILE_DIR, "examples", "basic"))
    calc.read_results()
    assert calc.results["forces"][0][0] == -1.3826330655069403


def test_qchem_read_results_intermediate():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    calc = QChem(atoms)
    os.chdir(os.path.join(FILE_DIR, "examples", "intermediate"))
    calc.read_results()
    assert calc.results["forces"][0][0] == -0.6955571014353796


def test_qchem_read_results_advanced():
    mol = Molecule.from_file(os.path.join(FILE_DIR, "test.xyz"))
    atoms = AseAtomsAdaptor.get_atoms(mol)
    calc = QChem(atoms)
    os.chdir(os.path.join(FILE_DIR, "examples", "advanced"))
    calc.read_results()
    assert calc.results["forces"][0][0] == -0.4270884974249971

