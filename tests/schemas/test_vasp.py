import os
from pathlib import Path
from shutil import which

import pytest
from ase.io import read
from monty.json import MontyDecoder, jsanitize

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "vasp_run1")


def mock_bader_analysis(*args, **kwargs):
    # NOTE: This is hard-coded for vasp_run1 with its 16 atoms
    bader_stats = {
        "min_dist": [1.0] * 16,
        "atomic_volume": [1.0] * 16,
        "vacuum_charge": 1.0,
        "vacuum_volume": 1.0,
        "bader_version": 1.0,
        "reference_used": [0.0] * 16,
        "partial_charges": [-1.0] * 16,
        "spin_moments": [0.0] * 16,
    }
    return bader_stats


@pytest.fixture(autouse=True)
def patch_pop_analyses(monkeypatch):
    # Monkeypatch the Bader analysis
    monkeypatch.setattr(
        "quacc.schemas.vasp.run_bader",
        mock_bader_analysis,
    )


def test_summarize_run():

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1)
    if results["nsites"] != len(atoms):
        raise AssertionError
    if results["atoms"] != atoms:
        raise AssertionError
    if results["output"]["energy"] != -33.15807349:
        raise AssertionError
    if results.get("calcs_reversed", None) is not None:
        raise AssertionError

    # Make sure default dir works
    cwd = os.getcwd()
    os.chdir(run1)
    summarize_run(atoms)
    os.chdir(cwd)

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1, compact=False, remove_empties=True)
    if results.get("calcs_reversed", None) is None:
        raise AssertionError
    if "author" in results:
        raise AssertionError
    if "additional_json" in results:
        raise AssertionError
    if "handler" in results["custodian"][0]:
        raise AssertionError
    if "corrections" in results["custodian"][0]:
        raise AssertionError

    # Make sure null are not removed
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1, remove_empties=False)
    if results["author"] != None:
        raise AssertionError
    if results["additional_json"] != {}:
        raise AssertionError
    if results["custodian"][0]["handler"] is not None:
        raise AssertionError
    if results["custodian"][0]["corrections"] != []:
        raise AssertionError

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms, dir_path=run1)
    results_atoms = results["atoms"]
    if atoms.info.get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError
    if results.get("atoms_info", {}) == {}:
        raise AssertionError
    if results["atoms_info"].get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError
    if results_atoms.info.get("test_dict", None) != {"hi": "there", "foo": "bar"}:
        raise AssertionError

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0, "magmoms": [2.0] * len(atoms)}
    results = summarize_run(atoms, dir_path=run1)
    results_atoms = results["atoms"]

    if atoms.calc is None:
        raise AssertionError
    if atoms.get_initial_magnetic_moments().tolist() != [3.14] * len(atoms):
        raise AssertionError

    if results_atoms.get_initial_magnetic_moments().tolist() != [2.0] * len(atoms):
        raise AssertionError
    if results_atoms.calc is not None:
        raise AssertionError

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, dir_path=run1, prep_next_run=False)
    if atoms.get_initial_magnetic_moments().tolist() != [3.14] * len(atoms):
        raise AssertionError
    results_atoms = results["atoms"]
    if results_atoms.get_initial_magnetic_moments().tolist() != [3.14] * len(atoms):
        raise AssertionError

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_bader_run():
    # Make sure Bader works
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, dir_path=run1, bader=True)
    struct = results["output"]["structure"]
    if struct.site_properties["bader_charge"] != [-1.0] * len(atoms):
        raise AssertionError
    if struct.site_properties["bader_spin"] != [0.0] * len(atoms):
        raise AssertionError
