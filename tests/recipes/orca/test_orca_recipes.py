from ase.build import molecule
from jobflow.managers.local import run_locally
import os
from pathlib import Path
from quacc.recipes.orca.core import RelaxMaker, StaticMaker
from shutil import copy

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
ORCA_DIR = os.path.join(FILE_DIR, "orca_run")


def setup_module():
    for f in os.listdir(ORCA_DIR):
        copy(os.path.join(ORCA_DIR, f), os.path.join(FILE_DIR, f))


def teardown_module():
    for f in os.listdir(ORCA_DIR):
        if os.path.exists(os.path.join(FILE_DIR, f)):
            os.remove(os.path.join(FILE_DIR, f))


def test_static_maker():

    atoms = molecule("H2")

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Static"

    job = StaticMaker(
        orcasimpleinput="HF def2-SVP SP", orcablock=r"%scf maxiter 300 end"
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Static"


def test_relax_maker():

    atoms = molecule("H2")

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Relax"

    job = RelaxMaker(
        orcasimpleinput="HF def2-SVP opt", orcablock=r"%scf maxiter 300 end"
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "ORCA-Relax"
