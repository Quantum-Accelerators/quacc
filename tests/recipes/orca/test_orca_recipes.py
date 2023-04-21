import multiprocessing
import os
from pathlib import Path
from shutil import copy, rmtree

from ase.build import molecule

from quacc.recipes.orca.core import RelaxJob, StaticJob

FILE_DIR = Path(__file__).resolve().parent
ORCA_DIR = os.path.join(FILE_DIR, "orca_run")


def setup_module():
    for f in os.listdir(ORCA_DIR):
        copy(os.path.join(ORCA_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(ORCA_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_Job():
    atoms = molecule("H2")
    nprocs = multiprocessing.cpu_count()

    output = StaticJob(atoms)
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp sp slowconv normalprint xyzfile"
    )
    assert output["parameters"]["orcablocks"] == f"%pal nprocs {nprocs} end"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1

    output = StaticJob(
        atoms,
        charge=-2,
        mult=3,
        input_swaps={"def2-SVP": True, "def2-TZVP": False},
        block_swaps={"%scf maxiter 300 end": True},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj sp slowconv normalprint xyzfile def2-svp"
    )
    assert (
        output["parameters"]["orcablocks"]
        == f"%scf maxiter 300 end %pal nprocs {nprocs} end"
    )


def test_relax_Job():
    atoms = molecule("H2")
    nprocs = multiprocessing.cpu_count()

    output = RelaxJob(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp opt slowconv normalprint xyzfile"
    )
    assert output["parameters"]["orcablocks"] == f"%pal nprocs {nprocs} end"

    output = RelaxJob(
        atoms,
        charge=-2,
        mult=3,
        input_swaps={
            "HF": True,
            "wb97x-d3bj": False,
            "def2-SVP": True,
            "def2-TZVP": False,
        },
        block_swaps={"%scf maxiter 300 end": True},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "opt slowconv normalprint xyzfile hf def2-svp"
    )
    assert (
        output["parameters"]["orcablocks"]
        == f"%scf maxiter 300 end %pal nprocs {nprocs} end"
    )
