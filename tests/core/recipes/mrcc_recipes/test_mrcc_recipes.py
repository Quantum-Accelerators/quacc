from __future__ import annotations

from pathlib import Path

import pytest
from ase import Atoms

from quacc.recipes.mrcc._base import prep_calculator, run_and_summarize
from quacc.recipes.mrcc.core import static_job

FILE_DIR = Path(__file__).parent


def test_static_job(tmp_path):
    atoms = Atoms("H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

    output = static_job(
        atoms,
        charge=0,
        spin_multiplicity=1,
        method="SCAN",
        basis="def2-SVP",
        mrccinput={"calc": "PBE", "basis": "STO-3G"},
        mrccblocks="symm=off",
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["mrccinput"] == {"basis": "STO-3G", "calc": "PBE"}
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0
    assert output["results"]["energy"] == pytest.approx(-2026.1497783941234)


def test_run_and_summarize(tmp_path):
    atoms = Atoms("H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

    output = run_and_summarize(
        atoms,
        charge=0,
        spin_multiplicity=1,
        default_inputs={"calc": "PBE", "basis": "STO-3G"},
        blocks="symm=off",
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["mrccinput"] == {"basis": "STO-3G", "calc": "PBE"}
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0
    assert output["results"]["energy"] == pytest.approx(-2026.1497783941234)


def test_prep_calculator():
    calc = prep_calculator(
        charge=2,
        spin_multiplicity=2,
        default_inputs={"calc": "HF", "basis": "def2-SVP"},
        blocks="symm=off",
        input_swaps={"calc": "PBE"},
    )

    ref_parameters = {
        "charge": 2,
        "mult": 2,
        "mrccinput": {"calc": "PBE", "basis": "def2-SVP"},
        "mrccblocks": "symm=off",
    }

    assert calc.parameters == ref_parameters
