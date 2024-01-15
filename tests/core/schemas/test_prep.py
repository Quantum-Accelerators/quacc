from copy import deepcopy
from pathlib import Path

import pytest
from ase import Atoms
from ase.build import bulk, molecule
from ase.io import read

from quacc.atoms.core import get_atoms_id
from quacc.calculators.vasp import Vasp
from quacc.schemas.prep import prep_magmoms, prep_next_run


@pytest.fixture()
def atoms_mag():
    file_dir = Path(__file__).parent
    return read(file_dir / ".." / "calculators" / "vasp" / "OUTCAR_mag.gz")


@pytest.fixture()
def atoms_nomag():
    file_dir = Path(__file__).parent
    return read(file_dir / ".." / "calculators" / "vasp" / "OUTCAR_nomag.gz")


@pytest.fixture()
def atoms_nospin():
    file_dir = Path(__file__).parent
    return read(file_dir / ".." / "calculators" / "vasp" / "OUTCAR_nospin.gz")


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


def test_prep_next_run(
    atoms_mag, atoms_nomag, atoms_nospin
):  # sourcery skip: extract-duplicate-method
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
    assert atoms.info.get("_old_ids", None) == [md5hash, md5hash]
    assert atoms.info.get("_id", None) == new_md5hash


def test_prep_magmoms():
    atoms = deepcopy(atoms_mag)
    atoms.info["test"] = "hi"
    mag = atoms.get_magnetic_moment()
    init_mags = atoms.get_initial_magnetic_moments()
    mags = atoms.get_magnetic_moments()
    atoms = prep_magmoms(atoms)
    assert atoms.info["test"] == "hi"
    assert atoms.calc is None
    assert atoms.get_initial_magnetic_moments().tolist() == mags.tolist()

    atoms = deepcopy(atoms_mag)
    atoms.info["test"] = "hi"
    atoms = prep_magmoms(atoms)
    assert atoms.info.get("test", None) == "hi"
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}

    atoms = deepcopy(atoms_nomag)
    mag = atoms.get_magnetic_moment()
    atoms = prep_magmoms(atoms)
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_magmoms(atoms)

    atoms = deepcopy(atoms_nospin)
    atoms = prep_magmoms(atoms)
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
