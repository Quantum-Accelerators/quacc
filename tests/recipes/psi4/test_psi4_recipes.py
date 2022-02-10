import os
from pathlib import Path
from shutil import copy

from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.gaussian.core import StaticMaker

FILE_DIR = Path(__file__).resolve().parent
PSI4_DIR = os.path.join(FILE_DIR, "psi4_run")


def setup_module():
    for f in os.listdir(PSI4_DIR):
        copy(os.path.join(PSI4_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(PSI4_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))


def test_static_maker():

    atoms = molecule("H2")

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Psi4-Static"
    assert "charge" not in output["parameters"]
    assert "mult" not in output["parameters"]
    assert output["parameters"]["method"] == "wb97x-d"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["num_thread"] == "max"

    job = StaticMaker(
        xc="m06l",
        basis="def2-svp",
        pop="regular",
        molden=False,
        swaps={"num_thread": 1, "mem": None},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert output["name"] == "Psi4-Static"
    assert output["parameters"]["method"] == "m06l"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["num_thread"] == 1
    assert "mem" not in output["parameters"]
