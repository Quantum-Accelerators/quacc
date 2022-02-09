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
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Gaussian-Static"
    assert output["parameters"]["sp"] == ""
    assert output["parameters"]["xc"] == "wB97X-D"
    assert output["parameters"]["basis"] == "def2-TZVP"
    assert output["parameters"]["integral"] == "ultrafine"
    assert output["parameters"]["gfinput"] == ""
    assert output["parameters"]["ioplist"] == ["6/7=3"]

    job = StaticMaker(
        xc="M06L",
        basis="def2-SVP",
        pop="regular",
        molden=False,
        swaps={"integral": "superfinegrid"},
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Gaussian-Static"
    assert output["parameters"]["sp"] == ""
    assert output["parameters"]["xc"] == "M06L"
    assert output["parameters"]["basis"] == "def2-SVP"
    assert output["parameters"]["integral"] == "superfinegrid"
    assert output["parameters"].get("gfinput", None) is None
    assert output["parameters"].get("ioplist", None) is None


def test_relax_maker():

    atoms = molecule("H2")

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Gaussian-Relax"
    assert output["parameters"]["opt"] == ""
    assert output["parameters"]["xc"] == "wB97X-D"
    assert output["parameters"]["basis"] == "def2-TZVP"
    assert output["parameters"]["integral"] == "ultrafine"

    job = RelaxMaker(
        xc="M06L",
        basis="def2-SVP",
        freq=True,
        swaps={"integral": "superfinegrid"},
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "Gaussian-Relax"
    assert output["parameters"]["opt"] == ""
    assert output["parameters"]["freq"] == ""
    assert output["parameters"]["xc"] == "M06L"
    assert output["parameters"]["basis"] == "def2-SVP"
    assert output["parameters"]["integral"] == "superfinegrid"
