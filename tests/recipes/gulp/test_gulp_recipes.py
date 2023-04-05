import os
from pathlib import Path
from shutil import copy, rmtree

from ase.build import molecule, bulk
from jobflow.managers.local import run_locally

from quacc.recipes.gulp.core import RelaxJob, StaticJob

FILE_DIR = Path(__file__).resolve().parent
GULP_DIR = os.path.join(FILE_DIR, "gulp_run")


def setup_module():
    for f in os.listdir(GULP_DIR):
        copy(os.path.join(GULP_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(GULP_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_Job():
    atoms = molecule("H2O")

    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "GULP-Static"

    atoms = bulk("Cu") * (2, 2, 2)

    job = StaticJob(gfnff=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "GULP-Static"


def test_relax_Job():
    atoms = molecule("H2O")

    job = RelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "GULP-Relax"

    atoms = bulk("Cu") * (2, 2, 2)

    job = RelaxJob(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "GULP-Relax"
