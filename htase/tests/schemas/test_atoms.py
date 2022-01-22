from htase.schemas.atoms import atoms_to_db
from ase.io.jsonio import decode
from ase.build import bulk


def test_atoms_to_db():

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_db(atoms)
    results_atoms = decode(results["atoms"])
    assert results_atoms.info.get("test", None) == "hi"

    atoms = bulk("Cu")
    atoms.info["test"] = "hi"
    results = atoms_to_db(atoms, strip_info=True)
    results_atoms = decode(results["atoms"])
    assert results_atoms == atoms
    assert results_atoms.info == {}
    assert results["nsites"] == len(atoms)
    assert results["atoms_info"].get("test", None) == "hi"

    atoms = bulk("Cu")
    results = atoms_to_db(atoms, get_metadata=False)
    assert decode(results["atoms"]) == atoms
    assert results.get("nsites", None) is None
