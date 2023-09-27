import os
from copy import deepcopy
from pathlib import Path

from ase import Atoms
from ase.build import bulk, molecule
from ase.io import read

from quacc.calculators.vasp import Vasp
from quacc.runners.prep import prep_next_run
from quacc.utils.atoms import get_atoms_id

FILE_DIR = Path(__file__).resolve().parent
ATOMS_MAG = read(os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nomag.gz")
)
ATOMS_NOSPIN = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nospin.gz")
)


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
    atoms = prep_next_run(atoms)
    assert atoms.info.get("test", None) == "hi"
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}

    atoms = deepcopy(ATOMS_MAG)
    atoms = prep_next_run(atoms, move_magmoms=False)
    assert atoms.get_initial_magnetic_moments().tolist() == init_mags.tolist()

    atoms = deepcopy(ATOMS_NOMAG)
    mag = atoms.get_magnetic_moment()
    atoms = prep_next_run(atoms)
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms)

    atoms = deepcopy(ATOMS_NOSPIN)
    atoms = prep_next_run(atoms)
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"magmom": mag - 2}
