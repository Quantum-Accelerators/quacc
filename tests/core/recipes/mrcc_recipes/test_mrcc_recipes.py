from __future__ import annotations

from pathlib import Path

import pytest
from ase import Atoms

from quacc import change_settings
from quacc.recipes.mrcc._base import prep_calculator, run_and_summarize
from quacc.recipes.mrcc.core import static_job

FILE_DIR = Path(__file__).parent


def test_static_job(tmp_path):
    atoms = Atoms("H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    with change_settings({"RESULTS_DIR": tmp_path}):
        output = static_job(
            atoms,
            charge=0,
            spin_multiplicity=1,
            method="SCAN",
            basis="def2-TZVP",
            calc="PBE",
        )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["basis"] == "def2-TZVP"
    assert output["parameters"]["calc"] == "PBE"
    assert output["parameters"]["symm"] == "off"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0
    assert output["results"]["energy"] == pytest.approx(-2061.4010013440234)

    # Check if it runs without specifying anything besides atoms
    atoms = Atoms("H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    with change_settings({"RESULTS_DIR": tmp_path}):
        output = static_job(atoms)

    assert output["results"]["energy"] == pytest.approx(-2061.4010013440234)


def test_run_and_summarize(tmp_path):
    atoms = Atoms("H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

    with change_settings({"RESULTS_DIR": tmp_path}):
        output = run_and_summarize(
            atoms,
            charge=0,
            spin_multiplicity=1,
            default_inputs={"calc": "PBE", "basis": "def2-tzvp", "symm": "off"},
        )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["calc"] == "PBE"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0
    assert output["results"]["energy"] == pytest.approx(-2061.4010013440234)


def test_prep_calculator():
    with pytest.raises(
        ValueError,
        match="For spin_multiplicity > 1, scftype keyword must be specified in mrccinput",
    ):
        calc = prep_calculator(
            spin_multiplicity=2,
            default_inputs={"calc": "PBE", "basis": "def2-tzvp", "symm": "off"},
        )

    with pytest.raises(
        ValueError,
        match="For spin_multiplicity > 1, scftype must not be set to RHF or RKS",
    ):
        calc = prep_calculator(
            spin_multiplicity=2,
            default_inputs={
                "calc": "HF",
                "basis": "def2-tzvp",
                "symm": "off",
                "scftype": "RHF",
            },
        )

    calc = prep_calculator(
        charge=2,
        spin_multiplicity=1,
        default_inputs={"calc": "HF", "basis": "def2-tzvp"},
        input_swaps={
            "calc": "PBE",
            "basis": "def2-SVP",
            "dfbasis_scf": """atomtype
H: def2-SVP
O: def2-SVP""",
        },
    )
    ref_parameters = {
        "calc": "PBE",
        "basis": "def2-SVP",
        "dfbasis_scf": "atomtype\nH: def2-SVP\nO: def2-SVP",
        "charge": 2,
        "mult": 1,
    }

    assert calc.parameters == ref_parameters
