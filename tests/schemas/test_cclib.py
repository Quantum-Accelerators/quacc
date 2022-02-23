import os
from pathlib import Path

import pytest
from ase.build import bulk
from ase.io import read
from monty.json import MontyDecoder, jsanitize

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.cclib import summarize_run

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "gaussian_run1")
log1 = os.path.join(run1, "gautest.log")


def test_summarize_run():

    # Make sure metadata is made
    atoms = read(log1)
    results = summarize_run(atoms, ".log", dir_path=run1)
    if results["nsites"] != len(atoms):
        raise AssertionError
    if results["atoms"] != atoms:
        raise AssertionError
    if results["spin_multiplicity"] != 1:
        raise AssertionError
    if results["nsites"] != 6:
        raise AssertionError
    if results["metadata"].get("success", None) != True:
        raise AssertionError

    # Make sure default dir works
    cwd = os.getcwd()
    os.chdir(run1)
    summarize_run(atoms, ".log")
    os.chdir(cwd)

    # Make sure info tags are handled appropriately
    atoms = read(log1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms, ".log", dir_path=run1)
    if atoms.info.get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError
    if results.get("atoms_info", {}) == {}:
        raise AssertionError
    if results["atoms_info"].get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError
    if results["atoms"].info.get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms, ".log", dir_path=run1)

    if atoms.calc is None:
        raise AssertionError
    if atoms.get_initial_magnetic_moments().tolist() != [3.14] * len(atoms):
        raise AssertionError

    if results["atoms"].get_initial_magnetic_moments().tolist() != [2.0] * len(
        atoms
    ):
        raise AssertionError
    if results["atoms"].calc is not None:
        raise AssertionError

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, ".log", dir_path=run1, prep_next_run=False)
    if atoms.get_initial_magnetic_moments().tolist() != [3.14] * len(atoms):
        raise AssertionError
    if results["atoms"].get_initial_magnetic_moments().tolist() != [3.14] * len(
        atoms
    ):
        raise AssertionError

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_errors():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        summarize_run(atoms, ".log", dir_path=run1)

    atoms = SmartVasp(atoms)
    with pytest.raises(ValueError):
        summarize_run(atoms, ".log", dir_path=run1)
