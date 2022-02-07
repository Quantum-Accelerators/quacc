import os
from pathlib import Path

from ase.io import read
from monty.json import MontyDecoder, jsanitize

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "vasp_run1")


def test_summarize_run():

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["output"]["energy"] == -33.15807349
    assert results.get("calcs_reversed", None) is None

    # Make sure default dir works
    cwd = os.getcwd()
    os.chdir(run1)
    summarize_run(atoms)
    os.chdir(cwd)

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1, compact=False, remove_empties=True)
    assert results.get("calcs_reversed", None) is not None
    assert "author" not in results
    assert "additional_json" not in results
    assert "handler" not in results["custodian"][0]
    assert "corrections" not in results["custodian"][0]

    # Make sure null are not removed
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1, remove_empties=False)
    assert results["author"] == None
    assert results["additional_json"] == {}
    assert results["custodian"][0]["handler"] is None
    assert results["custodian"][0]["corrections"] == []

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms, dir_path=run1)
    results_atoms = results["atoms"]
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results_atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0, "magmoms": [2.0] * len(atoms)}
    results = summarize_run(atoms, dir_path=run1)
    results_atoms = results["atoms"]

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results_atoms.get_initial_magnetic_moments().tolist() == [2.0] * len(atoms)
    assert results_atoms.calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, dir_path=run1, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    results_atoms = results["atoms"]
    assert results_atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
