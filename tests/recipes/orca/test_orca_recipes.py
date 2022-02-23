import multiprocessing
import os
from pathlib import Path
from shutil import copy

from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.orca.core import RelaxMaker, StaticMaker

FILE_DIR = Path(__file__).resolve().parent
ORCA_DIR = os.path.join(FILE_DIR, "orca_run")


def setup_module():
    for f in os.listdir(ORCA_DIR):
        copy(os.path.join(ORCA_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(ORCA_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))


def test_static_maker():

    atoms = molecule("H2")
    nprocs = multiprocessing.cpu_count()

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "ORCA-Static":
        raise AssertionError
    if (
        output["parameters"]["orcasimpleinput"] != "wb97x-d def2-tzvp sp slowconv normalprint"
    ):
        raise AssertionError
    if output["parameters"]["orcablocks"] != f"%pal nprocs {nprocs} end":
        raise AssertionError
    if output["parameters"]["charge"] != 0:
        raise AssertionError
    if output["parameters"]["mult"] != 1:
        raise AssertionError

    job = StaticMaker(
        input_swaps={"def2-SVP": True, "def2-TZVP": False},
        block_swaps={"%scf maxiter 300 end": True},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "ORCA-Static":
        raise AssertionError
    if output["parameters"]["charge"] != -2:
        raise AssertionError
    if output["parameters"]["mult"] != 3:
        raise AssertionError
    if (
        output["parameters"]["orcasimpleinput"] != "wb97x-d sp slowconv normalprint def2-svp"
    ):
        raise AssertionError
    if (
        output["parameters"]["orcablocks"] != f"%scf maxiter 300 end %pal nprocs {nprocs} end"
    ):
        raise AssertionError


def test_relax_maker():

    atoms = molecule("H2")
    nprocs = multiprocessing.cpu_count()

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["charge"] != 0:
        raise AssertionError
    if output["parameters"]["mult"] != 1:
        raise AssertionError
    if output["name"] != "ORCA-Relax":
        raise AssertionError
    if (
        output["parameters"]["orcasimpleinput"] != "wb97x-d def2-tzvp opt slowconv normalprint"
    ):
        raise AssertionError
    if output["parameters"]["orcablocks"] != f"%pal nprocs {nprocs} end":
        raise AssertionError

    job = RelaxMaker(
        input_swaps={
            "HF": True,
            "wb97x-d": False,
            "def2-SVP": True,
            "def2-TZVP": False,
        },
        block_swaps={"%scf maxiter 300 end": True},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "ORCA-Relax":
        raise AssertionError
    if output["parameters"]["charge"] != -2:
        raise AssertionError
    if output["parameters"]["mult"] != 3:
        raise AssertionError
    if (
        output["parameters"]["orcasimpleinput"] != "opt slowconv normalprint hf def2-svp"
    ):
        raise AssertionError
    if (
        output["parameters"]["orcablocks"] != f"%scf maxiter 300 end %pal nprocs {nprocs} end"
    ):
        raise AssertionError
