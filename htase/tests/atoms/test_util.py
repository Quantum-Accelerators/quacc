from ase.io import read
from ase.build import bulk
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from htase.util.atoms import (
    invert_atoms,
    make_slabs_from_bulk,
    make_max_slabs_from_bulk,
)
from htase.util.calc import cache_calc
from htase.calculators.vasp import SmartVasp
from ase.io.jsonio import encode, decode
from pathlib import Path
import os
import numpy as np
import pytest
from copy import deepcopy

FILE_DIR = Path(__file__).resolve().parent
ATOMS_MAG = read(os.path.join(FILE_DIR, "..", "vasp", "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(os.path.join(FILE_DIR, "..", "vasp", "OUTCAR_nomag.gz"))
ATOMS_NOSPIN = read(os.path.join(FILE_DIR, "..", "vasp", "OUTCAR_nospin.gz"))


def test_cache_calc():
    atoms = deepcopy(ATOMS_MAG)
    mag = atoms.get_magnetic_moment()
    atoms = cache_calc(atoms)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"magmom": mag - 2}
    atoms = cache_calc(atoms)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc1", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    assert atoms.info["results"]["calc1"]["magmom"] == mag - 2
    assert decode(encode(atoms)) == atoms

    atoms = deepcopy(ATOMS_NOMAG)
    mag = atoms.get_magnetic_moment()
    atoms = cache_calc(atoms)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"magmom": mag - 2}
    atoms = cache_calc(atoms)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc1", None) is not None
    assert atoms.info["results"]["calc0"]["magmom"] == mag
    assert atoms.info["results"]["calc1"]["magmom"] == mag - 2
    assert decode(encode(atoms)) == atoms

    atoms = deepcopy(ATOMS_NOSPIN)
    atoms = cache_calc(atoms)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"].get("magmom", None) is None
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"magmom": mag - 2}
    atoms = cache_calc(atoms)
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc1", None) is not None
    assert atoms.info["results"]["calc0"].get("magmom", None) is None
    assert atoms.info["results"]["calc1"]["magmom"] == mag - 2
    assert decode(encode(atoms)) == atoms


# def test_invert_atoms():
#     struct = Structure.from_file(os.path.join(FILE_DIR, "MnO2_primitive.cif.gz"))
#     slab = SlabGenerator(struct, [0, 0, 1], 10, 10).get_slab()
#     inverted_slab = invert_atoms(slab, return_struct=True)
#     assert pytest.approx(slab[0].x, 1e-9) == inverted_slab[0].x
#     assert pytest.approx(slab[0].y, 1e-9) == inverted_slab[0].y
#     assert pytest.approx(slab[0].z, 1e-5) == inverted_slab[0].z - 7.38042756

#     atoms = bulk("Cu")
#     struct = AseAtomsAdaptor.get_structure(atoms)
#     slab = [slab_struct for slab_struct in generate_all_slabs(struct, 1, 12, 20)][0]
#     true_slab = Structure.from_file(os.path.join(FILE_DIR, "slab_invert1.cif.gz"))
#     assert np.allclose(slab.frac_coords, true_slab.frac_coords)

#     inverted_slab = invert_atoms(slab, return_struct=True)
#     true_inverted_slab = Structure.from_file(
#         os.path.join(FILE_DIR, "slab_invert2.cif.gz")
#     )
#     assert np.allclose(inverted_slab.frac_coords, true_inverted_slab.frac_coords)


# This needs more tests
def test_make_slabs_from_bulk():

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_slabs_from_bulk(atoms)
    assert len(slabs) == 7
    shifts = [np.round(slab.info["slab_stats"]["shift"], 3) for slab in slabs]
    shifts.sort()
    assert shifts == [-0.875, -0.375, -0.125, 0.125, 0.25, 0.375, 0.875]

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_slabs_from_bulk(atoms, flip_asymmetric=False)
    assert len(slabs) == 4

    atoms = bulk("Cu")
    atoms.info["user_comments"] = "hi"
    slabs = make_slabs_from_bulk(atoms)
    assert slabs[-1].info.get("user_comments", None) == "hi"

    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms, required_surface_atoms=["Co"])
    assert slabs is None

    slabs = make_slabs_from_bulk(atoms, max_index=2)
    assert len(slabs) == 9

    slabs = make_slabs_from_bulk(atoms, z_fix=0.0)
    for slab in slabs:
        assert len(slab.constraints) == 0

    slabs = make_slabs_from_bulk(atoms, min_length_width=20)

    atoms = deepcopy(ATOMS_MAG)
    slabs = make_slabs_from_bulk(atoms)
    assert slabs[0].get_magnetic_moments()[0] == atoms.get_magnetic_moments()[0]
    assert slabs[-1].info.get("slab_stats", None) is not None


def make_max_slabs_from_bulk():
    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms)
    slabs2 = make_max_slabs_from_bulk(atoms, None)
    assert slabs == slabs2

    atoms = bulk("Cu")
    slabs = make_max_slabs_from_bulk(atoms, 2)
    assert len(slabs) == 2
    assert slabs[-1].info.get("slab_stats", None) is not None
