from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.json import unjsonify
from ase.build import bulk


def test_atoms_to_metadata():

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms)
    results_atoms = unjsonify(results["atoms"])
    assert results_atoms.info.get("test", None) == "hi"

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_metadata(atoms, strip_info=True)
    results_atoms = unjsonify(results["atoms"])
    assert results_atoms == atoms
    assert results_atoms.info == {}
    assert results["nsites"] == len(atoms)
    assert results["atoms_info"].get("test", None) == "hi"

    atoms = bulk("Cu")
    results = atoms_to_metadata(atoms, get_metadata=False)
    assert unjsonify(results["atoms"]) == atoms
    assert results.get("nsites", None) is None

    atoms = bulk("Cu")
    parent = bulk("Al") * (2, 1, 1)
    atoms.info = {"parent": parent}
    results = atoms_to_metadata(atoms)
    assert unjsonify(results["atoms_info"]["parent"]["atoms"]) == parent
    assert results["atoms_info"]["parent"]["nsites"] == len(parent)
