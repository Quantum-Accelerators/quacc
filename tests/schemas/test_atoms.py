from ase.build import bulk
from monty.json import MontyDecoder, jsanitize

from quacc.schemas.atoms import atoms_to_metadata


def test_atoms_to_metadata():

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    if results["atoms"].info.get("test", None) != "hi":
        raise AssertionError

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms, strip_info=True)
    if results["atoms"] != atoms:
        raise AssertionError
    if results["atoms"].info != {}:
        raise AssertionError
    if results["nsites"] != len(atoms):
        raise AssertionError
    if results["atoms_info"].get("test", None) != "hi":
        raise AssertionError

    atoms = bulk("Cu")
    results = atoms_to_metadata(atoms, get_metadata=False)
    if results["atoms"] != atoms:
        raise AssertionError
    if results.get("nsites", None) is not None:
        raise AssertionError

    atoms = bulk("Cu")
    parent = bulk("Al") * (2, 1, 1)
    atoms.info = {"parent": parent}
    results = atoms_to_metadata(atoms)
    if results["atoms_info"]["parent"]["atoms"] != parent:
        raise AssertionError
    if results["atoms_info"]["parent"]["nsites"] != len(parent):
        raise AssertionError

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
