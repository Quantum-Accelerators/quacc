import os

import pytest
from ase.build import bulk
from ase.optimize import BFGS
from jobflow.managers.local import run_locally

from quacc.recipes.emt.core import RelaxMaker, StaticMaker


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
        ):
            os.remove(f)


def test_static_maker():

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["asap_cutoff"] != False:
        raise AssertionError
    if output["name"] != "EMT-Static":
        raise AssertionError
    if output["results"]["energy"] != pytest.approx(0.07001766638245854):
        raise AssertionError

    job = StaticMaker(asap_cutoff=True).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["asap_cutoff"] != True:
        raise AssertionError
    if output["name"] != "EMT-Static":
        raise AssertionError
    if output["results"]["energy"] != pytest.approx(0.11074520235398744):
        raise AssertionError


def test_relax_maker():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["asap_cutoff"] != False:
        raise AssertionError
    if output["name"] != "EMT-Relax":
        raise AssertionError
    if output["results"]["energy"] != pytest.approx(-0.04517048198212592):
        raise AssertionError

    job = RelaxMaker(asap_cutoff=True).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["asap_cutoff"] != True:
        raise AssertionError
    if output["name"] != "EMT-Relax":
        raise AssertionError
    if output["results"]["energy"] != pytest.approx(-0.004527567070971017):
        raise AssertionError

    job = RelaxMaker(fmax=0.01).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "EMT-Relax":
        raise AssertionError
    if output["results"]["energy"] != pytest.approx(-0.0454470914411953):
        raise AssertionError

    job = RelaxMaker(optimizer=BFGS).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["name"] != "EMT-Relax":
        raise AssertionError
    if output["results"]["energy"] != pytest.approx(-0.0454470914411953):
        raise AssertionError
