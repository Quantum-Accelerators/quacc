from ase.io import read
from ase.build import bulk
from htase.calculators.vasp import SmartVasp
from pymatgen.core.surface import SlabGenerator
from pymatgen.core import Structure
from htase.util.atoms import (
    serialize,
    deserialize,
    make_conventional_cell,
    invert_slab,
    make_slabs_from_bulk,
)
from pathlib import Path
import os
import numpy as np
import pytest

FILE_DIR = Path(__file__).resolve().parent


def test_make_conventional_cell():
    atoms = read(os.path.join(FILE_DIR, "MnO2_primitive.cif.gz"))
    atoms = make_conventional_cell(atoms)
    truth = read(os.path.join(FILE_DIR, "MnO2_conventional.cif.gz"))
    assert np.array_equal(atoms.cell.lengths(), truth.cell.lengths())


def test_invert_slab():
    struct = Structure.from_file(os.path.join(FILE_DIR, "MnO2_primitive.cif.gz"))
    slab = SlabGenerator(struct, [0, 0, 1], 10, 10).get_slab()
    inverted_slab = invert_slab(slab, return_atoms=False)
    assert slab[0].x == inverted_slab[0].x
    assert slab[0].y == inverted_slab[0].y
    assert pytest.approx(slab[0].z, 1e-5) == inverted_slab[0].z - 7.38042756


def test_serialization():
    atoms = bulk("Cu")
    newatoms = deserialize(serialize(atoms))
    assert np.array_equal(atoms.get_positions(), newatoms.get_positions())

    atoms = read(os.path.join(FILE_DIR, "..", "vasp", "OUTCAR_mag.gz"))
    newatoms = deserialize(serialize(atoms))
    assert np.array_equal(atoms.get_positions(), newatoms.get_positions())
    assert atoms.get_potential_energy() == newatoms.get_potential_energy()
    assert np.array_equal(atoms.get_magnetic_moments(), newatoms.get_magnetic_moments())


# This test takes a while. Clearly, make_slabs_from_bulk could be improved
# for speed...
# def test_make_slabs_from_bulk():
#     atoms = bulk("Cu") * (3, 3, 3)
#     slabs = make_slabs_from_bulk(atoms)

#     slabs = make_slabs_from_bulk(atoms, max_index=0)

#     slabs = make_slabs_from_bulk(atoms, z_fix=0.0)

#     slabs = make_slabs_from_bulk(atoms, min_length_width=20)
