import os
from copy import deepcopy
from pathlib import Path

from ase.build import bulk, molecule
from ase.io import read

from quacc.calculators.vasp import SmartVasp
from quacc.util.atoms import (
    check_is_metal,
    get_atoms_id,
    get_highest_block,
    prep_next_run,
)

FILE_DIR = Path(__file__).resolve().parent
ATOMS_MAG = read(os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nomag.gz")
)
ATOMS_NOSPIN = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nospin.gz")
)


def test_get_atoms_id():
    atoms = bulk("Cu")
    md5hash = "d4859270a1a67083343bec0ab783f774"
    if get_atoms_id(atoms) != md5hash:
        raise AssertionError

    atoms.info["test"] = "hi"
    if get_atoms_id(atoms) != md5hash:
        raise AssertionError

    atoms.set_initial_magnetic_moments([1.0])
    md5maghash = "7d456a48c235e05cf17da4abcc433a4f"
    if get_atoms_id(atoms) != md5maghash:
        raise AssertionError


def test_prep_next_run():
    atoms = bulk("Cu")
    md5hash = "d4859270a1a67083343bec0ab783f774"
    atoms = prep_next_run(atoms)
    if atoms.info.get("_id", None) != md5hash:
        raise AssertionError
    if atoms.info.get("_old_ids", None) is not None:
        raise AssertionError
    atoms = prep_next_run(atoms)
    if atoms.info.get("_id", None) != md5hash:
        raise AssertionError
    if atoms.info.get("_old_ids", None) != [md5hash]:
        raise AssertionError
    atoms[0].symbol = "Pt"
    new_md5hash = "52087d50a909572d58e01cfb49d4911b"
    atoms = prep_next_run(atoms)
    if atoms.info.get("_old_ids", None) != [
        md5hash,
        md5hash,
    ]:
        raise AssertionError
    if atoms.info.get("_id", None) != new_md5hash:
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    atoms.info["test"] = "hi"
    mag = atoms.get_magnetic_moment()
    init_mags = atoms.get_initial_magnetic_moments()
    mags = atoms.get_magnetic_moments()
    atoms = prep_next_run(atoms)
    if atoms.info["test"] != "hi":
        raise AssertionError
    if atoms.calc != None:
        raise AssertionError
    if atoms.get_initial_magnetic_moments().tolist() != mags.tolist():
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    atoms.info["test"] = "hi"
    atoms = prep_next_run(atoms, store_results=True)
    if atoms.info.get("test", None) != "hi":
        raise AssertionError
    if atoms.info.get("results", None) is None:
        raise AssertionError
    if atoms.info["results"].get("calc0", None) is None:
        raise AssertionError
    if atoms.info["results"]["calc0"]["magmom"] != mag:
        raise AssertionError
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms, store_results=True)
    if atoms.info.get("results", None) is None:
        raise AssertionError
    if atoms.info["results"].get("calc1", None) is None:
        raise AssertionError
    if atoms.info["results"]["calc0"]["magmom"] != mag:
        raise AssertionError
    if atoms.info["results"]["calc1"]["magmom"] != mag - 2:
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    atoms = prep_next_run(atoms, move_magmoms=False)
    if atoms.get_initial_magnetic_moments().tolist() != init_mags.tolist():
        raise AssertionError

    atoms = deepcopy(ATOMS_NOMAG)
    mag = atoms.get_magnetic_moment()
    atoms = prep_next_run(atoms, store_results=True)
    if atoms.info.get("results", None) is None:
        raise AssertionError
    if atoms.info["results"].get("calc0", None) is None:
        raise AssertionError
    if atoms.info["results"]["calc0"]["magmom"] != mag:
        raise AssertionError
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms, store_results=True)
    if atoms.info.get("results", None) is None:
        raise AssertionError
    if atoms.info["results"].get("calc1", None) is None:
        raise AssertionError
    if atoms.info["results"]["calc0"]["magmom"] != mag:
        raise AssertionError
    if atoms.info["results"]["calc1"]["magmom"] != mag - 2:
        raise AssertionError

    atoms = deepcopy(ATOMS_NOSPIN)
    atoms = prep_next_run(atoms, store_results=True)
    if atoms.info.get("results", None) is None:
        raise AssertionError
    if atoms.info["results"].get("calc0", None) is None:
        raise AssertionError
    if atoms.info["results"]["calc0"].get("magmom", None) is not None:
        raise AssertionError
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"magmom": mag - 2}
    atoms = prep_next_run(atoms, store_results=True)
    if atoms.info.get("results", None) is None:
        raise AssertionError
    if atoms.info["results"].get("calc1", None) is None:
        raise AssertionError
    if atoms.info["results"]["calc0"].get("magmom", None) is not None:
        raise AssertionError
    if atoms.info["results"]["calc1"]["magmom"] != mag - 2:
        raise AssertionError


def test_check_is_metal():
    atoms = bulk("Cu")
    if check_is_metal(atoms) != True:
        raise AssertionError
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "O"
    if check_is_metal(atoms) != False:
        raise AssertionError
    atoms = molecule("H2O")
    if check_is_metal(atoms) != False:
        raise AssertionError


def test_get_highest_block():
    atoms = bulk("Cu")
    if get_highest_block(atoms) != "d":
        raise AssertionError
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "U"
    if get_highest_block(atoms) != "f":
        raise AssertionError
    atoms = molecule("H2O")
    if get_highest_block(atoms) != "p":
        raise AssertionError
