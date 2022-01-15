from ase.io import read
from ase.build import bulk, fcc100
from ase.io.jsonio import encode, decode
from ase.build import molecule
from htase.util.atoms import (
    flip_atoms,
    make_slabs_from_bulk,
    make_max_slabs_from_bulk,
    make_adsorbate_structures,
)
from htase.util.calc import cache_calc
from htase.calculators.vasp import SmartVasp
from pathlib import Path
import os
import numpy as np
from copy import deepcopy

FILE_DIR = Path(__file__).resolve().parent
ATOMS_MAG = read(os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nomag.gz")
)
ATOMS_NOSPIN = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nospin.gz")
)


def test_cache_calc():
    atoms = deepcopy(ATOMS_MAG)
    atoms.info = {"test": "hi"}
    mag = atoms.get_magnetic_moment()
    init_mags = atoms.get_initial_magnetic_moments()
    mags = atoms.get_magnetic_moments()
    atoms = cache_calc(atoms)
    assert atoms.info.get("test", None) == "hi"
    assert atoms.get_initial_magnetic_moments().tolist() == mags.tolist()
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

    atoms = deepcopy(ATOMS_MAG)
    atoms = cache_calc(atoms, move_magmoms=False)
    assert atoms.get_initial_magnetic_moments().tolist() == init_mags.tolist()

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


def test_flip_atoms():
    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    atoms.info["test"] = "hi"
    atoms.set_initial_magnetic_moments(
        [2.0 if atom.symbol == "Zn" else 1.0 for atom in atoms]
    )
    new_atoms = flip_atoms(atoms)
    assert (
        np.unique(np.round(atoms.get_all_distances(mic=True), 4)).tolist()
        == np.unique(np.round(new_atoms.get_all_distances(mic=True), 4)).tolist()
    )
    Zn_idx = [atom.index for atom in new_atoms if atom.symbol == "Zn"]
    Te_idx = [atom.index for atom in new_atoms if atom.symbol == "Te"]
    assert np.all(new_atoms.get_initial_magnetic_moments()[Zn_idx] == 2.0)
    assert np.all(new_atoms.get_initial_magnetic_moments()[Te_idx] == 1.0)
    assert new_atoms.info.get("test", None) == "hi"


def test_make_slabs_from_bulk():
    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    atoms.info["test"] = "hi"
    slabs = make_slabs_from_bulk(atoms)
    assert len(slabs) == 7
    assert len(slabs[0].constraints) != 0
    shifts = [np.round(slab.info["slab_stats"]["shift"], 3) for slab in slabs]
    shifts.sort()
    assert shifts == [-0.875, -0.375, -0.125, 0.125, 0.25, 0.375, 0.875]
    z_store = -np.inf
    for atom in slabs[3]:
        if atom.z > z_store:
            z_store = atom.z
            highest_atom = atom.symbol
    z_store = -np.inf
    for atom in slabs[6]:
        if atom.z > z_store:
            z_store = atom.z
            highest_atom2 = atom.symbol
    assert highest_atom != highest_atom2
    assert atoms.info.get("test", None) == "hi"

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_slabs_from_bulk(atoms, flip_asymmetric=False)
    assert len(slabs) == 4

    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms, allowed_surface_atoms=["Co"])
    assert slabs is None

    slabs = make_slabs_from_bulk(atoms, max_index=2)
    assert len(slabs) == 9

    slabs = make_slabs_from_bulk(atoms, z_fix=0.0)
    for slab in slabs:
        assert len(slab.constraints) == 0

    slabs = make_slabs_from_bulk(atoms, min_length_width=20)
    for slab in slabs:
        assert slab.cell.lengths()[0] >= 20
        assert slab.cell.lengths()[1] >= 20

    atoms = deepcopy(ATOMS_MAG)
    slabs = make_slabs_from_bulk(atoms)
    assert slabs[0].get_magnetic_moments()[0] == atoms.get_magnetic_moments()[0]
    assert slabs[-1].info.get("slab_stats", None) is not None

    atoms = read(os.path.join(FILE_DIR, "Zn2CuAu.cif.gz"))
    min_d = atoms.get_all_distances(mic=True)
    min_d = np.min(min_d[min_d != 0.0])
    slabs = make_slabs_from_bulk(atoms)
    assert len(slabs) == 31
    assert (
        slabs[3].info["slab_stats"]["shift"] == -slabs[22].info["slab_stats"]["shift"]
    )
    z_store = -np.inf
    for atom in slabs[3]:
        if atom.z > z_store:
            z_store = atom.z
            highest_atom = atom.symbol
    z_store = -np.inf
    for atom in slabs[22]:
        if atom.z > z_store:
            z_store = atom.z
            highest_atom2 = atom.symbol
    assert highest_atom != highest_atom2
    # This is to make sure nothing funky happened to our atom positions...
    for slab in slabs:
        d = slab.get_all_distances(mic=True)
        assert np.round(np.min(d[d != 0]), 4) == np.round(min_d, 4)


def test_make_max_slabs_from_bulk():
    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms)
    slabs2 = make_max_slabs_from_bulk(atoms, None)
    assert slabs == slabs2

    atoms = bulk("Cu")
    slabs = make_max_slabs_from_bulk(atoms, 2)
    assert len(slabs) == 2
    assert slabs[-1].info.get("slab_stats", None) is not None

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_max_slabs_from_bulk(atoms, 4)
    assert len(slabs) == 4


def test_make_adsorbate_structures():

    atoms = fcc100("Cu", size=(2, 2, 2))
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)
    new_atoms = make_adsorbate_structures(atoms, "H2O", modes=["ontop"])
    assert new_atoms[0].has("initial_magmoms") is False

    mol = molecule("O2")
    mol.set_initial_magnetic_moments([1.0, 1.0])
    new_atoms = make_adsorbate_structures(atoms, mol)
    assert len(new_atoms) == 3
    assert new_atoms[0].get_initial_magnetic_moments().tolist() == [0.0] * len(
        atoms
    ) + [1.0, 1.0]

    atoms = fcc100("Cu", size=(2, 2, 2))
    mags = [5.0] * len(atoms)
    atoms.set_initial_magnetic_moments(mags)
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)
    mol = molecule("O2")
    mol.set_initial_magnetic_moments([1.0, 1.0])
    new_atoms = make_adsorbate_structures(atoms, mol)
    assert len(new_atoms) == 3
    assert new_atoms[0].get_initial_magnetic_moments().tolist() == mags + [1.0, 1.0]

    new_atoms = make_adsorbate_structures(atoms, "H2O")
    assert len(new_atoms) == 3
    assert new_atoms[0].get_initial_magnetic_moments().tolist() == mags + [0, 0, 0]
    new_atoms = make_adsorbate_structures(atoms, "H2O", modes=["ontop"])
    assert len(new_atoms) == 1

    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_symbols=["Cu", "Fe"]
    )
    assert len(new_atoms) == 3

    new_atoms = make_adsorbate_structures(atoms, "H2O", allowed_surface_indices=[23])
    assert len(new_atoms) == 2

    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_indices=[18], modes=["ontop"]
    )
    assert len(new_atoms) == 1

    atoms[18].symbol = "Fe"
    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_symbols=["Fe"], modes=["ontop"]
    )
    assert len(new_atoms) == 1
    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_symbols=["Cu", "Fe"]
    )
    assert len(new_atoms) == 10

    assert new_atoms[0].info.get("adsorbates", None) is not None
    assert decode(new_atoms[0].info.get("adsorbates", None)[0]["atoms"]) == molecule(
        "H2O"
    )
