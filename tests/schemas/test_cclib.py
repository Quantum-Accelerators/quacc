import os
from pathlib import Path

import pytest
from ase.build import bulk
from ase.io import read
from monty.json import MontyDecoder, jsanitize

from quacc.calculators.vasp import Vasp
from quacc.schemas.cclib import summarize_run

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "gaussian_run1")
log1 = os.path.join(run1, "Gaussian.log")


def test_summarize_run():
    # Make sure metadata is made
    atoms = read(log1)
    results = summarize_run(atoms, ".log", dir_path=run1, check_convergence=True)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["spin_multiplicity"] == 1
    assert results["natoms"] == 6
    assert results["metadata"].get("success", None) == True
    assert "pull_request" in results["builder_meta"]

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)

    # Make sure default dir works
    cwd = os.getcwd()
    os.chdir(run1)
    summarize_run(atoms, ".log")
    os.chdir(cwd)

    # Make sure info tags are handled appropriately
    atoms = read(log1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms, ".log", dir_path=run1)
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms, ".log", dir_path=run1)

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results["atoms"].get_initial_magnetic_moments().tolist() == [2.0] * len(
        atoms
    )
    assert results["atoms"].calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, ".log", dir_path=run1, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    assert results["atoms"].get_initial_magnetic_moments().tolist() == [3.14] * len(
        atoms
    )

    # Test remove_empties
    # Make sure metadata is made
    atoms = read(log1)
    results = summarize_run(atoms, ".log", dir_path=run1, remove_empties=True)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["spin_multiplicity"] == 1
    assert results["natoms"] == 6
    assert results["metadata"].get("success", None) == True
    assert "pull_request" not in results["builder_meta"]


def test_errors():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        summarize_run(atoms, ".log", dir_path=run1)

    calc = Vasp(atoms)
    atoms.calc = calc
    with pytest.raises(ValueError):
        summarize_run(atoms, ".log", dir_path=run1)
