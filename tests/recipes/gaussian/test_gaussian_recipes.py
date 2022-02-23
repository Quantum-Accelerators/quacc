import os
from pathlib import Path
from shutil import copy

from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.gaussian.core import RelaxMaker, StaticMaker

FILE_DIR = Path(__file__).resolve().parent
GAUSSIAN_DIR = os.path.join(FILE_DIR, "gaussian_run")


def setup_module():
    for f in os.listdir(GAUSSIAN_DIR):
        copy(os.path.join(GAUSSIAN_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(GAUSSIAN_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))


def test_static_maker():

    atoms = molecule("H2")

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "Gaussian-Static":
        raise AssertionError
    if "charge" in output["parameters"]:
        raise AssertionError
    if "mult" in output["parameters"]:
        raise AssertionError
    if output["parameters"]["sp"] != "":
        raise AssertionError
    if output["parameters"]["xc"] != "wb97x-d":
        raise AssertionError
    if output["parameters"]["basis"] != "def2-tzvp":
        raise AssertionError
    if output["parameters"]["integral"] != "ultrafine":
        raise AssertionError
    if output["parameters"]["gfinput"] != "":
        raise AssertionError
    if output["parameters"]["ioplist"] != ["6/7=3"]:
        raise AssertionError

    job = StaticMaker(
        xc="m06l",
        basis="def2-svp",
        pop="regular",
        molden=False,
        swaps={"integral": "superfinegrid"},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["charge"] != -2:
        raise AssertionError
    if output["parameters"]["mult"] != 3:
        raise AssertionError
    if output["name"] != "Gaussian-Static":
        raise AssertionError
    if output["parameters"]["sp"] != "":
        raise AssertionError
    if output["parameters"]["xc"] != "m06l":
        raise AssertionError
    if output["parameters"]["basis"] != "def2-svp":
        raise AssertionError
    if output["parameters"]["integral"] != "superfinegrid":
        raise AssertionError
    if "gfinput" in output["parameters"]:
        raise AssertionError
    if "ioptlist" in output["parameters"]:
        raise AssertionError
    if "opt" in output["parameters"]:
        raise AssertionError


def test_relax_maker():

    atoms = molecule("H2")

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "Gaussian-Relax":
        raise AssertionError
    if "charge" in output["parameters"]:
        raise AssertionError
    if "mult" in output["parameters"]:
        raise AssertionError
    if output["parameters"]["opt"] != "":
        raise AssertionError
    if output["parameters"]["xc"] != "wb97x-d":
        raise AssertionError
    if output["parameters"]["basis"] != "def2-tzvp":
        raise AssertionError
    if output["parameters"]["integral"] != "ultrafine":
        raise AssertionError
    if "freq" in output["parameters"]:
        raise AssertionError
    if "sp" in output["parameters"]:
        raise AssertionError

    job = RelaxMaker(
        xc="m06l", basis="def2-svp", freq=True, swaps={"integral": "superfinegrid"}
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["charge"] != -2:
        raise AssertionError
    if output["parameters"]["mult"] != 3:
        raise AssertionError
    if output["name"] != "Gaussian-Relax":
        raise AssertionError
    if output["parameters"]["opt"] != "":
        raise AssertionError
    if output["parameters"]["freq"] != "":
        raise AssertionError
    if output["parameters"]["xc"] != "m06l":
        raise AssertionError
    if output["parameters"]["basis"] != "def2-svp":
        raise AssertionError
    if output["parameters"]["integral"] != "superfinegrid":
        raise AssertionError
