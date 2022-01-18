import os
from htase.schemas.vasp.summarize import get_results
from htase.calculators.vasp import SmartVasp
from ase.io import read
from ase.io.jsonio import decode
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "run1")


def test_summarize():

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    results = get_results(dir_path=run1)
    assert results["nsites"] == len(atoms)

    # Make sure Atoms object is stored properly
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    results = get_results(dir_path=run1, atoms=atoms)
    assert decode(results["atoms"]) == atoms

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = get_results(dir_path=run1, atoms=atoms)
    results_atoms = decode(results["atoms"])
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results_atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0, "magmoms": [2.0] * len(atoms)}
    results = get_results(dir_path=run1, atoms=atoms)
    results_atoms = decode(results["atoms"])

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results_atoms.get_initial_magnetic_moments().tolist() == [2.0] * len(atoms)
    assert results_atoms.calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = get_results(dir_path=run1, atoms=atoms, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    results_atoms = decode(results["atoms"])
    assert results_atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
