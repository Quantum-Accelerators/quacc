import os
from copy import deepcopy
from pathlib import Path

from ase import Atoms
from ase.build import bulk, molecule
from ase.io import read
import pytest

from quacc.calculators.vasp import Vasp
from quacc.util.atoms import (
    check_is_metal,
    get_atoms_id,
    get_highest_block,
    prep_next_run,
    check_charge_and_spin,
)

FILE_DIR = Path(__file__).resolve().parent
ATOMS_MAG = read(os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nomag.gz")
)
ATOMS_NOSPIN = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nospin.gz")
)
OS_ATOMS = read(os.path.join(FILE_DIR, "OS_test.xyz"))


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


def test_prep_next_run():  # sourcery skip: extract-duplicate-method
    atoms = bulk("Cu")
    md5hash = "d4859270a1a67083343bec0ab783f774"
    atoms = prep_next_run(atoms)
    assert atoms.info.get("_id", None) == md5hash
    assert atoms.info.get("_old_ids", None) is None
    atoms = prep_next_run(atoms)
    assert atoms.info.get("_id", None) == md5hash
    assert atoms.info.get("_old_ids", None) == [md5hash]
    atoms[0].symbol = "Pt"
    new_md5hash = "52087d50a909572d58e01cfb49d4911b"
    atoms = prep_next_run(atoms)
    assert atoms.info.get("_old_ids", None) == [
        md5hash,
        md5hash,
    ]
    assert atoms.info.get("_id", None) == new_md5hash

    atoms = deepcopy(ATOMS_MAG)
    atoms.info["test"] = "hi"
    mag = atoms.get_magnetic_moment()
    init_mags = atoms.get_initial_magnetic_moments()
    mags = atoms.get_magnetic_moments()
    atoms = prep_next_run(atoms)
    assert atoms.info["test"] == "hi"
    assert atoms.calc is None
    assert atoms.get_initial_magnetic_moments().tolist() == mags.tolist()

    atoms = deepcopy(ATOMS_MAG)
    atoms.info["test"] = "hi"
    atoms = prep_next_run(atoms, store_results=True)
    assert atoms.info.get("test", None) == "hi"
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms, store_results=True)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc1", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    assert atoms.info["results"]["calc1"]["magmom"] == mag - 2

    atoms = deepcopy(ATOMS_MAG)
    atoms = prep_next_run(atoms, move_magmoms=False)
    assert atoms.get_initial_magnetic_moments().tolist() == init_mags.tolist()

    atoms = deepcopy(ATOMS_NOMAG)
    mag = atoms.get_magnetic_moment()
    atoms = prep_next_run(atoms, store_results=True)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms, store_results=True)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc1", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    assert atoms.info["results"]["calc1"]["magmom"] == mag - 2

    atoms = deepcopy(ATOMS_NOSPIN)
    atoms = prep_next_run(atoms, store_results=True)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"].get("magmom", None) is None
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms, store_results=True)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc1", None) is not None
    assert atoms.info["results"]["calc0"].get("magmom", None) is None
    assert atoms.info["results"]["calc1"]["magmom"] == mag - 2


def test_check_is_metal():
    atoms = bulk("Cu")
    assert check_is_metal(atoms) is True
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "O"
    assert check_is_metal(atoms) is False
    atoms = molecule("H2O")
    assert check_is_metal(atoms) is False


def test_get_highest_block():
    atoms = bulk("Cu")
    assert get_highest_block(atoms) == "d"
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "U"
    assert get_highest_block(atoms) == "f"
    atoms = molecule("H2O")
    assert get_highest_block(atoms) == "p"


def test_check_charge_and_spin():
    atoms = molecule("CH3")
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == 0
    assert spin_multiplicity == 2
    with pytest.raises(RuntimeError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, spin_multiplicity=1)
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, charge=0, spin_multiplicity=1)
    with pytest.raises(RuntimeError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, spin_multiplicity=3)
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(atoms, charge=0, spin_multiplicity=3)
    charge, spin_multiplicity = check_charge_and_spin(atoms, charge=-1)
    assert charge == -1
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(atoms, charge=-1, spin_multiplicity=3)
    assert charge == -1
    assert spin_multiplicity == 3
    charge, spin_multiplicity = check_charge_and_spin(OS_ATOMS)
    assert charge == 0
    assert spin_multiplicity == 2
    charge, spin_multiplicity = check_charge_and_spin(OS_ATOMS, charge=1)
    assert charge == 1
    assert spin_multiplicity == 1
    charge, spin_multiplicity = check_charge_and_spin(OS_ATOMS, charge=0, spin_multiplicity=4)
    assert charge == 0
    assert spin_multiplicity == 4
    with pytest.raises(ValueError):
        charge, spin_multiplicity = check_charge_and_spin(OS_ATOMS, charge=0, spin_multiplicity=3)
    atoms = molecule("O2")
    charge, spin_multiplicity = check_charge_and_spin(atoms)
    assert charge == 0
    assert spin_multiplicity == 3
    charge, spin_multiplicity = check_charge_and_spin(atoms, charge=0)
    assert charge == 0
    assert spin_multiplicity == 1 # Somewhat controvertial if this should end up at 1 or 3
