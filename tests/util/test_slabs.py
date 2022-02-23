import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, fcc100, molecule
from ase.io import read

from quacc.util.slabs import (
    flip_atoms,
    make_adsorbate_structures,
    make_max_slabs_from_bulk,
    make_slabs_from_bulk,
)

FILE_DIR = Path(__file__).resolve().parent
ATOMS_MAG = read(os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nomag.gz")
)
ATOMS_NOSPIN = read(
    os.path.join(FILE_DIR, "..", "calculators", "vasp", "OUTCAR_nospin.gz")
)


def test_flip_atoms():
    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    atoms.info["test"] = "hi"
    atoms.set_initial_magnetic_moments(
        [2.0 if atom.symbol == "Zn" else 1.0 for atom in atoms]
    )
    new_atoms = flip_atoms(atoms)
    if (
        np.unique(np.round(atoms.get_all_distances(mic=True), 4)).tolist() != np.unique(np.round(new_atoms.get_all_distances(mic=True), 4)).tolist()
    ):
        raise AssertionError
    Zn_idx = [atom.index for atom in new_atoms if atom.symbol == "Zn"]
    Te_idx = [atom.index for atom in new_atoms if atom.symbol == "Te"]
    if not np.all(new_atoms.get_initial_magnetic_moments()[Zn_idx] == 2.0):
        raise AssertionError
    if not np.all(new_atoms.get_initial_magnetic_moments()[Te_idx] == 1.0):
        raise AssertionError
    if new_atoms.info.get("test", None) != "hi":
        raise AssertionError


def test_make_slabs_from_bulk():
    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    atoms.info["test"] = "hi"
    slabs = make_slabs_from_bulk(atoms)
    if len(slabs) != 7:
        raise AssertionError
    if len(slabs[0].constraints) == 0:
        raise AssertionError
    shifts = [np.round(slab.info["slab_stats"]["shift"], 3) for slab in slabs]
    shifts.sort()
    if shifts != [-0.875, -0.375, -0.125, 0.125, 0.25, 0.375, 0.875]:
        raise AssertionError
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
    if highest_atom == highest_atom2:
        raise AssertionError
    if atoms.info.get("test", None) != "hi":
        raise AssertionError

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_slabs_from_bulk(atoms, flip_asymmetric=False)
    if len(slabs) != 4:
        raise AssertionError

    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms, allowed_surface_atoms=["Co"])
    if slabs is not None:
        raise AssertionError

    slabs = make_slabs_from_bulk(atoms, max_index=2)
    if len(slabs) != 9:
        raise AssertionError

    slabs = make_slabs_from_bulk(atoms, z_fix=0.0)
    for slab in slabs:
        if len(slab.constraints) != 0:
            raise AssertionError

    slabs = make_slabs_from_bulk(atoms, min_length_width=20)
    for slab in slabs:
        if slab.cell.lengths()[0] < 20:
            raise AssertionError
        if slab.cell.lengths()[1] < 20:
            raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    slabs = make_slabs_from_bulk(atoms)
    if slabs[0].get_magnetic_moments()[0] != atoms.get_magnetic_moments()[0]:
        raise AssertionError
    if slabs[-1].info.get("slab_stats", None) is None:
        raise AssertionError

    atoms = read(os.path.join(FILE_DIR, "Zn2CuAu.cif.gz"))
    min_d = atoms.get_all_distances(mic=True)
    min_d = np.min(min_d[min_d != 0.0])
    slabs = make_slabs_from_bulk(atoms)
    if len(slabs) != 31:
        raise AssertionError
    if (
        slabs[3].info["slab_stats"]["shift"] != -slabs[22].info["slab_stats"]["shift"]
    ):
        raise AssertionError
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
    if highest_atom == highest_atom2:
        raise AssertionError
    # This is to make sure nothing funky happened to our atom positions...
    for slab in slabs:
        d = slab.get_all_distances(mic=True)
        if np.round(np.min(d[d != 0]), 4) != np.round(min_d, 4):
            raise AssertionError


def test_make_max_slabs_from_bulk():
    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms)
    slabs2 = make_max_slabs_from_bulk(atoms, None)
    if slabs != slabs2:
        raise AssertionError

    atoms = bulk("Cu")
    slabs = make_max_slabs_from_bulk(atoms, 2)
    if len(slabs) != 2:
        raise AssertionError
    if slabs[-1].info.get("slab_stats", None) is None:
        raise AssertionError

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_max_slabs_from_bulk(atoms, 4)
    if len(slabs) != 4:
        raise AssertionError


def test_make_adsorbate_structures():

    atoms = fcc100("Cu", size=(2, 2, 2))
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)
    new_atoms = make_adsorbate_structures(atoms, "H2O", modes=["ontop"])
    if new_atoms[0].has("initial_magmoms") is not False:
        raise AssertionError

    mol = molecule("O2")
    mol.set_initial_magnetic_moments([1.0, 1.0])
    new_atoms = make_adsorbate_structures(atoms, mol)
    if len(new_atoms) != 3:
        raise AssertionError
    if new_atoms[0].get_initial_magnetic_moments().tolist() != [0.0] * len(
        atoms
    ) + [1.0, 1.0]:
        raise AssertionError

    atoms = fcc100("Cu", size=(2, 2, 2))
    mags = [5.0] * len(atoms)
    atoms.set_initial_magnetic_moments(mags)
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)
    mol = molecule("O2")
    mol.set_initial_magnetic_moments([1.0, 1.0])
    new_atoms = make_adsorbate_structures(atoms, mol)
    if len(new_atoms) != 3:
        raise AssertionError
    if new_atoms[0].get_initial_magnetic_moments().tolist() != mags + [1.0, 1.0]:
        raise AssertionError

    new_atoms = make_adsorbate_structures(atoms, "H2O")
    if len(new_atoms) != 3:
        raise AssertionError
    if new_atoms[0].get_initial_magnetic_moments().tolist() != mags + [0, 0, 0]:
        raise AssertionError
    new_atoms = make_adsorbate_structures(atoms, "H2O", modes=["ontop"])
    if len(new_atoms) != 1:
        raise AssertionError

    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_symbols=["Cu", "Fe"]
    )
    if len(new_atoms) != 3:
        raise AssertionError

    new_atoms = make_adsorbate_structures(atoms, "H2O", allowed_surface_indices=[6])
    if len(new_atoms) != 2:
        raise AssertionError

    atoms[7].symbol = "Fe"
    new_atoms = make_adsorbate_structures(atoms, "H2O", modes=["ontop"])
    if len(new_atoms) != 3:
        raise AssertionError

    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_symbols=["Fe"], modes=["ontop"]
    )
    if len(new_atoms) != 1:
        raise AssertionError

    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_indices=[7], modes=["ontop"]
    )
    if len(new_atoms) != 1:
        raise AssertionError

    new_atoms = make_adsorbate_structures(
        atoms, "H2O", allowed_surface_symbols=["Cu", "Fe"]
    )
    if len(new_atoms) != 6:
        raise AssertionError
    if new_atoms[0].info.get("adsorbates", None) is None:
        raise AssertionError
    if new_atoms[0].info["adsorbates"][0]["adsorbate"] != molecule("H2O"):
        raise AssertionError


def test_errors():
    atoms = fcc100("Cu", size=(2, 2, 2))
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)

    with pytest.raises(ValueError):
        make_adsorbate_structures(
            atoms, "H2O", min_distance=1.0, find_ads_sites_kwargs={"distance": 1.0}
        )
    with pytest.raises(ValueError):
        make_adsorbate_structures(
            atoms,
            "H2O",
            modes=["ontop"],
            find_ads_sites_kwargs={"positions": ["ontop"]},
        )
    with pytest.raises(ValueError):
        make_adsorbate_structures(atoms, "H2O", allowed_surface_indices=[100])
    with pytest.raises(ValueError):
        make_adsorbate_structures(atoms, "WOW")
