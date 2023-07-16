import os
from pathlib import Path
from shutil import copy

from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.gaussian.core import relax_job, static_job

FILE_DIR = Path(__file__).resolve().parent
GAUSSIAN_DIR = os.path.join(FILE_DIR, "gaussian_run")


DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    SETTINGS.CREATE_UNIQUE_WORKDIR = False


def teardown_module():
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR


def prep_files():
    for f in os.listdir(GAUSSIAN_DIR):
        copy(os.path.join(GAUSSIAN_DIR, f), f)


def test_static_job(tmpdir):
    tmpdir.chdir()
    prep_files()

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


def test_relax_job(tmpdir):
    tmpdir.chdir()
    prep_files()

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
