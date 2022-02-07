from ase.build import bulk
from monty.json import MontyDecoder, jsanitize

from quacc.schemas.atoms import atoms_to_metadata


def test_atoms_to_metadata():

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    assert results["atoms"].info.get("test", None) == "hi"

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms, strip_info=True)
    assert results["atoms"] == atoms
    assert results["atoms"].info == {}
    assert results["nsites"] == len(atoms)
    assert results["atoms_info"].get("test", None) == "hi"

    atoms = bulk("Cu")
    results = atoms_to_metadata(atoms, get_metadata=False)
    assert results["atoms"] == atoms
    assert results.get("nsites", None) is None

    atoms = bulk("Cu")
    parent = bulk("Al") * (2, 1, 1)
    atoms.info = {"parent": parent}
    results = atoms_to_metadata(atoms)
    assert results["atoms_info"]["parent"]["atoms"] == parent
    assert results["atoms_info"]["parent"]["nsites"] == len(parent)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
