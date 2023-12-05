from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.io import read
from monty.json import MontyDecoder, jsanitize
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.schemas.atoms import atoms_to_metadata


@pytest.fixture()
def test_cifs():
    file_dir = Path(__file__).parent

    return file_dir / "test_files"


def test_atoms_to_metadata(test_cifs):
    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    assert results["atoms"].info.get("test", None) == "hi"
    assert results["structure"] == AseAtomsAdaptor.get_structure(atoms)
    assert "molecule" not in results
    assert "pymatgen_version" in results["builder_meta"]

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "X"
    atoms.info["test"] = "hi"
    atoms_no_dummy = atoms.copy()
    del atoms_no_dummy[[atom.index for atom in atoms if atom.symbol == "X"]]
    results = atoms_to_metadata(atoms)
    assert results["atoms"].info.get("test", None) == "hi"
    assert results["structure"] == AseAtomsAdaptor.get_structure(atoms_no_dummy)

    atoms = bulk("Cu")
    results = atoms_to_metadata(atoms, get_metadata=False)
    assert results["atoms"] == atoms
    assert results.get("nsites") is None

    atoms = molecule("H2O")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    assert "structure" not in results
    assert results["atoms"].info.get("test", None) == "hi"
    assert results["molecule"] == AseAtomsAdaptor().get_molecule(atoms)

    atoms = bulk("Cu")
    parent = bulk("Al") * (2, 1, 1)
    atoms.info = {
        "parent": parent,
        "test": [
            parent,
            1.0,
            [2.0, 3.0],
            {"subtest": parent, "subtest2": "hi"},
            np.array([1.0, 2.0]),
        ],
    }
    results = atoms_to_metadata(atoms)
    assert results["atoms_info"]["parent"]["atoms"] == parent
    assert results["atoms_info"]["parent"]["nsites"] == len(parent)
    assert results["atoms_info"]["test"][0]["atoms"] == parent
    assert results["atoms_info"]["test"][0]["nsites"] == len(parent)
    assert results["atoms_info"]["test"][1] == 1.0
    assert results["atoms_info"]["test"][2] == [2.0, 3.0]
    assert results["atoms_info"]["test"][3]["subtest"]["atoms"] == parent
    assert results["atoms_info"]["test"][3]["subtest"]["nsites"] == len(parent)
    assert results["atoms_info"]["test"][3]["subtest2"] == "hi"
    assert results["atoms_info"]["test"][4] == [1.0, 2.0]

    atoms = read(test_cifs / "nonserializable_info.cif.gz")
    atoms.info["parent"] = parent
    results = atoms_to_metadata(atoms)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
