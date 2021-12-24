import os
from ase.io import read
from ase.io.jsonio import encode, decode
from htase.calculators.vasp import SmartVasp
from pathlib import Path
import numpy as np

FILE_DIR = Path(__file__).resolve().parent


def test_magmom_carryover():
    atoms = read(os.path.join(FILE_DIR, "OUTCAR_mag.gz"))
    assert atoms.get_magnetic_moments()[0] == 0.468
    atoms = decode(encode(atoms))
    assert atoms.get_magnetic_moments()[0] == 0.468
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    assert atoms.get_initial_magnetic_moments()[0] == 0.468

    atoms = read(os.path.join(FILE_DIR, "OUTCAR_nomag.gz"))
    assert atoms.get_magnetic_moments()[0] == 0.468
    atoms = decode(encode(atoms))
    assert atoms.get_magnetic_moments()[0] == 0.468
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    assert atoms.has("initial_magmoms") is False
    assert np.all(atoms.get_initial_magnetic_moments() == 0)
