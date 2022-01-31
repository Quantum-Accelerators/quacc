import os
from ase.io import read
from quacc.schemas.calc import summarize_run
from quacc.util.json import unjsonify
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "vasp_run1")


def test_summarize_run():

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms)
    assert results["nsites"] == len(atoms)
    assert unjsonify(results["atoms"]) == atoms

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms)
    results_atoms = unjsonify(results["atoms"])
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results_atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms)
    results_atoms = unjsonify(results["atoms"])

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results_atoms.get_initial_magnetic_moments().tolist() == [2.0] * len(atoms)
    assert results_atoms.calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    results_atoms = unjsonify(results["atoms"])
    assert results_atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
