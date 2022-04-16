import os

import pytest
from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

from quacc.recipes.dftb.core import StaticJob


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or f.endswith(".bin")
            or f.endswith(".hsd")
            or f.endswith(".gz")
            or f.endswith(".gen")
            or f.endswith(".tag")
        ):
            os.remove(f)


def test_static_Job():

    atoms = molecule("H2O")

    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Static"
    assert output["results"]["energy"] == pytest.approx(-137.97459675495645)

    atoms = bulk("Cu")

    job = StaticJob(kpts=(3, 3, 3)).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Static"
    assert output["results"]["energy"] == pytest.approx(-107.55154244254307)
