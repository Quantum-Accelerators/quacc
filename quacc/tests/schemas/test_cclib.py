import os
from quacc.schemas.cclib import summarize_run
from ase.io import read
from ase.io.jsonio import decode
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "gaussian_run1")
log1 = os.path.join(run1, "gautest.log.gz")


def test_summarize_run():

    # Make sure metadata is made
    atoms = read(log1)
    results = summarize_run(atoms, ".log", dir_path=FILE_DIR)
    assert results["nsites"] == len(atoms)
    assert decode(results["atoms"]) == atoms
    assert results["attributes"].get("mult", None) == 1
    assert results["attributes"].get("natom", None) == 6
    assert results["metadata"].get("success", None) == True

    # Make sure info tags are handled appropriately
    atoms = read(log1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms, ".log", dir_path=FILE_DIR)
    results_atoms = decode(results["atoms"])
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results_atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms, ".log", dir_path=FILE_DIR)
    results_atoms = decode(results["atoms"])

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results_atoms.get_initial_magnetic_moments().tolist() == [2.0] * len(atoms)
    assert results_atoms.calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, ".log", dir_path=FILE_DIR, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    results_atoms = decode(results["atoms"])
    assert results_atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
