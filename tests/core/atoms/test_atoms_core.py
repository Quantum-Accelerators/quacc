from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.atoms import Atoms
from ase.build import bulk, molecule
from ase.io import read

from quacc.atoms.core import check_charge_and_spin, check_is_metal, get_atoms_id

FILE_DIR = Path(__file__).parent


@pytest.fixture()
def atoms_mag():
    return read(FILE_DIR / ".." / "calculators" / "vasp" / "OUTCAR_mag.gz")


@pytest.fixture()
def atoms_nomag():
    return read(FILE_DIR / ".." / "calculators" / "vasp" / "OUTCAR_nomag.gz")


@pytest.fixture()
def atoms_nospin():
    return read(FILE_DIR / ".." / "calculators" / "vasp" / "OUTCAR_nospin.gz")


@pytest.fixture()
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


def test_check_is_metal():
    atoms = bulk("Cu")
    assert check_is_metal(atoms) is True
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "O"
    assert check_is_metal(atoms) is False
    atoms = molecule("H2O")
    assert check_is_metal(atoms) is False


def test_check_charge_and_spin(os_atoms):
    atoms = Atoms.fromdict(
        {
            "numbers": np.array([6, 1, 1, 1]),
            "positions": np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.0, 1.07841, 0.0],
                    [0.93393, -0.539205, 0.0],
                    [-0.93393, -0.539205, 0.0],
                ]
            ),
            "cell": np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            "pbc": np.array([False, False, False]),
        }
    )
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == 0
    assert spin_multiplicity == 2
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, spin_multiplicity=1)
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(
            atoms, charge=0, spin_multiplicity=1
        )
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, spin_multiplicity=3)
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(
            atoms, charge=0, spin_multiplicity=3
        )
    charge, spin_multiplicity = check_charge_and_spin(atoms, charge=-1)
    assert charge == -1
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(
        atoms, charge=-1, spin_multiplicity=3
    )
    assert charge == -1
    assert spin_multiplicity == 3
    charge, spin_multiplicity = check_charge_and_spin(os_atoms)
    assert charge == 0
    assert spin_multiplicity == 2
    charge, spin_multiplicity = check_charge_and_spin(os_atoms, charge=1)
    assert charge == 1
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(
        os_atoms, charge=0, spin_multiplicity=4
    )
    assert charge == 0
    assert spin_multiplicity == 4
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(
            os_atoms, charge=0, spin_multiplicity=3
        )

    atoms = molecule("CH3")
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == 0
    assert spin_multiplicity == 2
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, spin_multiplicity=1)
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(
            atoms, charge=0, spin_multiplicity=1
        )
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, spin_multiplicity=3)
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(
            atoms, charge=0, spin_multiplicity=3
        )
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, charge=-1)

    charge, spin_multiplicity = check_charge_and_spin(
        atoms, charge=-1, spin_multiplicity=3
    )
    assert charge == -1
    assert spin_multiplicity == 3
    charge, spin_multiplicity = check_charge_and_spin(os_atoms)
    assert charge == 0
    assert spin_multiplicity == 2
    charge, spin_multiplicity = check_charge_and_spin(os_atoms, charge=1)
    assert charge == 1
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(
        os_atoms, charge=0, spin_multiplicity=4
    )
    assert charge == 0
    assert spin_multiplicity == 4
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(
            os_atoms, charge=0, spin_multiplicity=3
        )

    atoms = molecule("CH3")
    atoms.charge = -2
    atoms.spin_multiplicity = 4
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == -2
    assert spin_multiplicity == 4

    atoms = molecule("CH3")
    atoms.set_initial_charges([-2, 0, 0, 0])
    atoms.set_initial_magnetic_moments([3, 0, 0, 0])
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == -2
    assert spin_multiplicity == 4

    atoms = molecule("CH3")
    atoms.set_initial_charges([-2, 0, 0, 0])
    atoms.set_initial_magnetic_moments([4, 0, 0, 0])
    with pytest.raises(ValueError):
        check_charge_and_spin(atoms)

    atoms = Atoms.fromdict(
        {
            "numbers": np.array([8, 8]),
            "positions": np.array([[0.0, 0.0, 0.622978], [0.0, 0.0, -0.622978]]),
            "cell": np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            "pbc": np.array([False, False, False]),
        }
    )
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == 0
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(atoms, charge=0)
    assert charge == 0
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(
        atoms, charge=0, spin_multiplicity=3
    )
    assert charge == 0
    assert spin_multiplicity == 3
    with pytest.raises(ValueError):
        check_charge_and_spin(atoms, charge=0, spin_multiplicity=2)

    atoms = molecule("O2")
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == 0
    assert spin_multiplicity == 3
    charge, spin_multiplicity = check_charge_and_spin(atoms, charge=0)
    assert charge == 0
    assert spin_multiplicity == 3
