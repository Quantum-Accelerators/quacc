import os
from pathlib import Path

import pytest
from ase.build import bulk
from ase.io import read
from monty.json import MontyDecoder, jsanitize

from quacc.schemas.calc import summarize_run

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "vasp_run1")


def test_summarize_run():

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms)
    if results["nsites"] != len(atoms):
        raise AssertionError
    if results["atoms"] != atoms:
        raise AssertionError

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms)
    if atoms.info.get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError
    if results.get("atoms_info", {}) == {}:
        raise AssertionError
    if results["atoms_info"].get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError
    if results["atoms"].info.get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms)

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
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, prep_next_run=False)
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
        summarize_run(atoms)

    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.calc.results = {}
    with pytest.raises(ValueError):
        summarize_run(atoms)
