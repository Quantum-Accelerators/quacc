import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, fcc100, molecule
from ase.io import read

from quacc.utils.slabs import (
    flip_atoms,
    get_surface_energy,
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
    slabs = make_slabs_from_bulk(atoms, allowed_surface_symbols=["Co"])
    assert slabs == []

    atoms = bulk("Cu")
    slabs = make_slabs_from_bulk(atoms, allowed_surface_symbols="Co")
    assert slabs == []

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
    slabs = make_slabs_from_bulk(atoms)
    slabs2 = make_max_slabs_from_bulk(atoms, None, allowed_surface_symbols="Cu")
    assert slabs == slabs2

    atoms = bulk("Cu")
    slabs = make_max_slabs_from_bulk(atoms, max_slabs=2)
    assert len(slabs) == 2
    assert slabs[-1].info.get("slab_stats", None) is not None

    atoms = bulk("Cu")
    slabs = make_max_slabs_from_bulk(atoms, max_slabs=2, randomize=True)
    assert len(slabs) == 2
    assert slabs[-1].info.get("slab_stats", None) is not None

    atoms = read(os.path.join(FILE_DIR, "ZnTe.cif.gz"))
    slabs = make_max_slabs_from_bulk(atoms, max_slabs=4)
    assert len(slabs) == 4


def test_make_adsorbate_structures():
    h2o = molecule("H2O")
    atoms = fcc100("Cu", size=(2, 2, 2))
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)
    new_atoms = make_adsorbate_structures(atoms, h2o, modes=["ontop"])
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

    new_atoms = make_adsorbate_structures(atoms, h2o)
    assert len(new_atoms) == 3
    assert new_atoms[0].get_initial_magnetic_moments().tolist() == mags + [0, 0, 0]
    new_atoms = make_adsorbate_structures(atoms, h2o, modes="ontop")
    assert len(new_atoms) == 1

    new_atoms = make_adsorbate_structures(
        atoms, h2o, allowed_surface_symbols=["Cu", "Fe"]
    )
    assert len(new_atoms) == 3

    new_atoms = make_adsorbate_structures(atoms, h2o, allowed_surface_indices=6)
    assert len(new_atoms) == 2

    atoms[7].symbol = "Fe"
    new_atoms = make_adsorbate_structures(atoms, h2o, modes=["ontop"])
    assert len(new_atoms) == 3

    new_atoms = make_adsorbate_structures(
        atoms, h2o, allowed_surface_symbols="Fe", modes=["ontop"]
    )
    assert len(new_atoms) == 1
    new_atoms = make_adsorbate_structures(
        atoms, h2o, allowed_surface_indices=[7], modes=["ontop"]
    )
    assert len(new_atoms) == 1

    new_atoms = make_adsorbate_structures(
        atoms, h2o, allowed_surface_symbols=["Cu", "Fe"]
    )
    assert len(new_atoms) == 6
    assert new_atoms[0].info.get("adsorbates", None) is not None
    assert new_atoms[0].info["adsorbates"][0]["adsorbate"] == h2o


def test_get_surface_energy():
    atoms = bulk("Cu")
    slab = fcc100("Cu", size=(2, 2, 2))
    bulk_energy = -1.0
    slab_energy = -20.0
    cleave_energy = get_surface_energy(atoms, slab, bulk_energy, slab_energy)
    assert cleave_energy == -0.23020081184152974


def test_errors():
    h2o = molecule("H2O")
    atoms = fcc100("Cu", size=(2, 2, 2))
    atoms.set_tags(None)
    atoms.center(vacuum=10, axis=2)

    with pytest.raises(ValueError):
        make_adsorbate_structures(
            atoms, h2o, min_distance=1.0, find_ads_sites_kwargs={"distance": 1.0}
        )
    with pytest.raises(ValueError):
        make_adsorbate_structures(
            atoms,
            h2o,
            modes=["ontop"],
            find_ads_sites_kwargs={"positions": ["ontop"]},
        )
    with pytest.raises(ValueError):
        make_adsorbate_structures(atoms, h2o, allowed_surface_indices=[100])
