import os
from shutil import rmtree, which

import numpy as np
import pytest
from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

from quacc.recipes.dftb.core import RelaxJob, StaticJob

DFTBPLUS_EXISTS = bool(which("dftb+"))


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
        if "quacc-tmp" in f or f == "tmp_dir":
            rmtree(f)


@pytest.mark.skipif(
    DFTBPLUS_EXISTS is False,
    reason="DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
def test_static_Job():

    atoms = molecule("H2O")

    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Static"
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.97459675495645)

    atoms = bulk("Cu")

    job = StaticJob(kpts=(3, 3, 3)).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Static"
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "3 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 3 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 3"
    assert output["results"]["energy"] == pytest.approx(-107.55154244254307)


@pytest.mark.skipif(
    DFTBPLUS_EXISTS is False,
    reason="DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
def test_relax_job():

    atoms = molecule("H2O")

    job = RelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Relax"
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.98336133558837)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is False
    )

    atoms = bulk("Cu")

    job = RelaxJob(kpts=(3, 3, 3)).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Relax"
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "3 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 3 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 3"
    assert output["parameters"]["Driver_"] == "GeometryOptimization"
    assert output["parameters"]["Driver_LatticeOpt"] == "No"
    assert output["atoms"].get_volume() == atoms.get_volume()

    atoms = bulk("Cu")

    job = RelaxJob(method="GFN1-xTB", kpts=(3, 3, 3), lattice_opt=True).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["name"] == "DFTB-Relax"
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN1-xTB"
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "3 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 3 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 3"
    assert output["parameters"]["Driver_"] == "GeometryOptimization"
    assert output["parameters"]["Driver_LatticeOpt"] == "Yes"
    assert output["atoms"].get_volume() != atoms.get_volume()
