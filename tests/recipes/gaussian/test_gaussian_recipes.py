import os
from pathlib import Path
from shutil import copy, rmtree

from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.gaussian.core import RelaxJob, StaticJob

FILE_DIR = Path(__file__).resolve().parent
GAUSSIAN_DIR = os.path.join(FILE_DIR, "gaussian_run")


def setup_module():
    for f in os.listdir(GAUSSIAN_DIR):
        copy(os.path.join(GAUSSIAN_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(GAUSSIAN_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f:
            rmtree(f)


def test_static_Job():

    atoms = molecule("H2")

    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Gaussian-Static"
    assert "charge" not in output["parameters"]
    assert "mult" not in output["parameters"]
    assert output["parameters"]["sp"] == ""
    assert output["parameters"]["xc"] == "wb97x-d"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["integral"] == "ultrafine"
    assert output["parameters"]["gfinput"] == ""
    assert output["parameters"]["ioplist"] == ["6/7=3"]

    job = StaticJob(
        xc="m06l",
        basis="def2-svp",
        pop="regular",
        molden=False,
        swaps={"integral": "superfinegrid"},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert output["name"] == "Gaussian-Static"
    assert output["parameters"]["sp"] == ""
    assert output["parameters"]["xc"] == "m06l"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["integral"] == "superfinegrid"
    assert "gfinput" not in output["parameters"]
    assert "ioptlist" not in output["parameters"]
    assert "opt" not in output["parameters"]


def test_relax_Job():

    atoms = molecule("H2")

    job = RelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Gaussian-Relax"
    assert "charge" not in output["parameters"]
    assert "mult" not in output["parameters"]
    assert output["parameters"]["opt"] == ""
    assert output["parameters"]["xc"] == "wb97x-d"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["integral"] == "ultrafine"
    assert "freq" not in output["parameters"]
    assert "sp" not in output["parameters"]

    job = RelaxJob(
        xc="m06l", basis="def2-svp", freq=True, swaps={"integral": "superfinegrid"}
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert output["name"] == "Gaussian-Relax"
    assert output["parameters"]["opt"] == ""
    assert output["parameters"]["freq"] == ""
    assert output["parameters"]["xc"] == "m06l"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["integral"] == "superfinegrid"
