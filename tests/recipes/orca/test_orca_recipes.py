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
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Static"
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d def2-tzvp sp slowconv normalprint"
    )
    assert output["parameters"]["orcablocks"] == f"%pal nprocs {nprocs} end"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1

    job = StaticMaker(
        input_swaps={"def2-SVP": True, "def2-TZVP": False},
        block_swaps={"%scf maxiter 300 end": True},
    ).make(atoms, charge=-2, mult=3)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Static"
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d sp slowconv normalprint def2-svp"
    )
    assert (
        output["parameters"]["orcablocks"] == r"%scf maxiter 300 end %pal nprocs 16 end"
    )


def test_relax_maker():

    atoms = molecule("H2")
    nprocs = multiprocessing.cpu_count()

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["name"] == "ORCA-Relax"
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d def2-tzvp opt slowconv normalprint"
    )
    assert output["parameters"]["orcablocks"] == f"%pal nprocs {nprocs} end"

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
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Relax"
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "opt slowconv normalprint hf def2-svp"
    )
    assert (
        output["parameters"]["orcablocks"] == r"%scf maxiter 300 end %pal nprocs 16 end"
    )
