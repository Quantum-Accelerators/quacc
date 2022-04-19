import os
from pathlib import Path

from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.psi4.core import StaticMaker

FILE_DIR = Path(__file__).resolve().parent


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if f.endswith(".dat"):
            os.remove(f)


def test_static_maker():

    atoms = molecule("H2")

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Psi4-Static"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["multiplicity"] == 1
    assert output["parameters"]["method"] == "wb97x-v"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["num_threads"] == "max"

    job = StaticMaker(
        method="m06l",
        basis="def2-svp",
        swaps={"num_threads": 1, "mem": None, "pop": "regular"},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["multiplicity"] == 3
    assert output["name"] == "Psi4-Static"
    assert output["parameters"]["method"] == "m06l"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["num_threads"] == 1
    assert output["parameters"]["pop"] == "regular"
    assert "mem" not in output["parameters"]
