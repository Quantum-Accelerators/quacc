from ase.io import read
from ase.build import bulk
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
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
    assert np.allclose(atoms.cell.lengths(), truth.cell.lengths())


def test_invert_slab():
    struct = Structure.from_file(os.path.join(FILE_DIR, "MnO2_primitive.cif.gz"))
    slab = SlabGenerator(struct, [0, 0, 1], 10, 10).get_slab()
    inverted_slab = invert_slab(slab, return_struct=True)
    assert slab[0].x == inverted_slab[0].x
    assert slab[0].y == inverted_slab[0].y
    assert pytest.approx(slab[0].z, 1e-5) == inverted_slab[0].z - 7.38042756

    atoms = bulk("Cu")
    struct = AseAtomsAdaptor.get_structure(atoms)
    slab = [slab_struct for slab_struct in generate_all_slabs(struct, 1, 12, 20)][0]
    true_slab = Structure.from_file(os.path.join(FILE_DIR, "slab_invert1.cif.gz"))
    assert np.allclose(slab.frac_coords, true_slab.frac_coords)

    inverted_slab = invert_slab(slab, return_struct=True)
    true_inverted_slab = Structure.from_file(
        os.path.join(FILE_DIR, "slab_invert2.cif.gz")
    )
    assert np.allclose(inverted_slab.frac_coords, true_inverted_slab.frac_coords)


def test_serialization():
    atoms = bulk("Cu")
    newatoms = deserialize(serialize(atoms))
    assert np.allclose(atoms.get_positions(), newatoms.get_positions())

    atoms = read(os.path.join(FILE_DIR, "..", "vasp", "OUTCAR_mag.gz"))
    newatoms = deserialize(serialize(atoms))
    assert np.allclose(atoms.get_positions(), newatoms.get_positions())
    assert atoms.get_potential_energy() == newatoms.get_potential_energy()
    assert np.allclose(atoms.get_magnetic_moments(), newatoms.get_magnetic_moments())

    atoms = read(os.path.join(FILE_DIR, "..", "vasp", "OUTCAR_nomag.gz"))
    newatoms = deserialize(serialize(atoms))
    assert np.allclose(atoms.get_magnetic_moments(), newatoms.get_magnetic_moments())


# This test takes a while. Clearly, make_slabs_from_bulk could be improved
# for speed...
# def test_make_slabs_from_bulk():
#     atoms = bulk("Cu")
#     slabs = make_slabs_from_bulk(atoms)

#     slabs = make_slabs_from_bulk(atoms, max_index=0)

#     slabs = make_slabs_from_bulk(atoms, z_fix=0.0)

#     slabs = make_slabs_from_bulk(atoms, min_length_width=20)
