from __future__ import annotations

from pathlib import Path

import pytest
from ase.build import bulk, molecule
from monty.json import MontyDecoder, jsanitize
from pymatgen.core import Molecule, Structure

from quacc.schemas.atoms import atoms_to_metadata


@pytest.fixture
def test_cifs():
    file_dir = Path(__file__).parent

    return file_dir / "test_files"


def test_atoms_to_metadata():
    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    assert results["atoms"].info.get("test", None) == "hi"
    assert results["structure"] == Structure.from_ase_atoms(atoms)
    assert "molecule" not in results
    assert "pymatgen_version" in results["builder_meta"]

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "X"
    atoms.info["test"] = "hi"
    atoms_no_dummy = atoms.copy()
    del atoms_no_dummy[[atom.index for atom in atoms if atom.symbol == "X"]]
    results = atoms_to_metadata(atoms)
    assert results["atoms"].info.get("test", None) == "hi"
    assert results["structure"] == Structure.from_ase_atoms(atoms_no_dummy)

    atoms = bulk("Cu")
    results = atoms_to_metadata(atoms, get_metadata=False)
    assert results["atoms"] == atoms
    assert results.get("nsites") is None

    atoms = molecule("H2O")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    assert "structure" not in results
    assert results["atoms"].info.get("test", None) == "hi"
    assert results["molecule"] == Molecule.from_ase_atoms(atoms)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
