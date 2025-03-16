from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase.atoms import Atoms
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.io import read
from numpy.testing import assert_allclose

from quacc.atoms.core import (
    check_is_metal,
    get_atoms_id,
    get_atoms_id_parsl,
    get_spin_multiplicity_attribute,
    perturb,
)

FILE_DIR = Path(__file__).parent


@pytest.fixture
def os_atoms():
    return read(FILE_DIR / "OS_test.xyz")


def test_init():
    atoms = bulk("Cu")
    assert Atoms.from_dict(atoms.as_dict()) == atoms

    atoms = molecule("CH3")
    assert Atoms.from_dict(atoms.as_dict()) == atoms


def test_get_atoms_id():
    atoms = bulk("Cu")
    md5hash = "d4859270a1a67083343bec0ab783f774"
    assert get_atoms_id(atoms) == md5hash

    atoms.info["test"] = "hi"
    assert get_atoms_id(atoms) == md5hash

    atoms.set_initial_magnetic_moments([1.0])
    md5maghash = "7d456a48c235e05cf17da4abcc433a4f"
    assert get_atoms_id(atoms) == md5maghash


def test_get_atoms_id_parsl():
    atoms = bulk("Cu")

    md5hash = b"\xd4\x85\x92p\xa1\xa6p\x834;\xec\n\xb7\x83\xf7t"
    assert get_atoms_id_parsl(atoms) == md5hash


def test_check_is_metal():
    atoms = bulk("Cu")
    assert check_is_metal(atoms) is True
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "O"
    assert check_is_metal(atoms) is False
    atoms = molecule("H2O")
    assert check_is_metal(atoms) is False


def test_perturb():
    atoms = Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)])
    atoms.get_positions().copy()
    matrix = [[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]
    scale = 0.5
    perturbed_atoms = perturb(atoms, matrix, scale)
    assert_allclose(
        perturbed_atoms.get_positions(), [[0.05, 0.05, 0.05], [0.1, 0.1, 0.84]]
    )
    assert atoms == Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)])


def test_get_spin_mult():
    atoms = molecule("H2O")
    assert get_spin_multiplicity_attribute(atoms) is None

    atoms = molecule("H2O")
    atoms.spin_multiplicity = 3
    assert get_spin_multiplicity_attribute(atoms) == 3

    atoms = molecule("H2O")
    atoms.set_initial_magnetic_moments([0.0, 0.0, 0.0])
    assert get_spin_multiplicity_attribute(atoms) == 1

    atoms = molecule("H2O")
    atoms.set_initial_magnetic_moments([0.0, 0.0, 0.0])
    atoms.calc = EMT()
    assert get_spin_multiplicity_attribute(atoms) == 1

    atoms = molecule("H2O")
    atoms.calc = EMT()
    atoms.get_potential_energy()
    atoms.calc.results["magmom"] = 1.999
    assert get_spin_multiplicity_attribute(atoms) == 3

    atoms = molecule("H2O")
    atoms.calc = EMT()
    atoms.get_potential_energy()
    atoms.calc.results["magmoms"] = np.array([0.0, 0.0, 1.999])
    assert get_spin_multiplicity_attribute(atoms) == 3
