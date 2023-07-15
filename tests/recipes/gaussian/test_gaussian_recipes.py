import os
from pathlib import Path
from shutil import copy, rmtree

from ase.build import molecule

from quacc.recipes.gaussian.core import relax_job, static_job

FILE_DIR = Path(__file__).resolve().parent
GAUSSIAN_DIR = os.path.join(FILE_DIR, "gaussian_run")


def setup_module():
    for f in os.listdir(GAUSSIAN_DIR):
        copy(os.path.join(GAUSSIAN_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(GAUSSIAN_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))


def teardown_function():
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_job():
    atoms = molecule("H2")

    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["sp"] == ""
    assert output["parameters"]["xc"] == "wb97x-d"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["integral"] == "ultrafine"
    assert output["parameters"]["gfinput"] == ""
    assert output["parameters"]["ioplist"] == [
        "6/7=3",
        "2/9=2000",
    ]  # see ASE issue #660

    output = static_job(
        atoms,
        charge=-2,
        multiplicity=3,
        xc="m06l",
        basis="def2-svp",
        calc_swaps={"integral": "superfinegrid"},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert output["parameters"]["sp"] == ""
    assert output["parameters"]["xc"] == "m06l"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["integral"] == "superfinegrid"
    assert output["parameters"]["gfinput"] == ""
    assert output["parameters"]["ioplist"] == [
        "6/7=3",
        "2/9=2000",
    ]  # see ASE issue #660
    assert "opt" not in output["parameters"]


def test_relax_Job():
    atoms = molecule("H2")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["opt"] == ""
    assert output["parameters"]["xc"] == "wb97x-d"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["integral"] == "ultrafine"
    assert "freq" not in output["parameters"]
    assert "sp" not in output["parameters"]

    output = relax_job(
        atoms,
        charge=-2,
        multiplicity=3,
        xc="m06l",
        basis="def2-svp",
        freq=True,
        calc_swaps={"integral": "superfinegrid"},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert output["parameters"]["opt"] == ""
    assert output["parameters"]["freq"] == ""
    assert output["parameters"]["xc"] == "m06l"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["integral"] == "superfinegrid"
    assert output["parameters"]["ioplist"] == ["2/9=2000"]  # see ASE issue #660
