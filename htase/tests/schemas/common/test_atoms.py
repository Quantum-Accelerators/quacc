from htase.schemas.common.atoms import atoms_to_db
from ase.io.jsonio import decode
from ase.build import bulk


def test_atoms_to_db():

    atoms = bulk("Cu")
    md5hash = "5f9a7d971f5b6f655ccbde7403f5b2fe"
    atoms.info["test"] = "hi"
    results = atoms_to_db(atoms)
    results_atoms = decode(results["atoms"])
    assert results_atoms.info.get("test", None) == "hi"
    assert results_atoms.info.get("_id", None) == md5hash
    assert results_atoms.info.get("_old_ids", None) is None
    results = atoms_to_db(results_atoms)
    results_atoms = decode(results["atoms"])
    assert results_atoms.info.get("_id", None) == md5hash
    assert results_atoms.info.get("_old_ids", None) == [md5hash]
    results_atoms[0].symbol = "Pt"
    new_md5hash = "31b258b71510a59c44ea4274eeb07a64"
    results = atoms_to_db(results_atoms)
    results_atoms = decode(results["atoms"])
    assert results_atoms.info.get("_old_ids", None) == [
        md5hash,
        md5hash,
    ]
    assert results_atoms.info.get("_id", None) == new_md5hash

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
