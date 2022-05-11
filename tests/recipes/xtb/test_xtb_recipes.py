import os

import numpy as np
import pytest
from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

try:
    from xtb.ase.calculator import XTB
except (ModuleNotFoundError, ImportError):
    XTB = None
from quacc.recipes.xtb.core import RelaxJob, StaticJob, ThermoJob


def teardown_module():
    for f in os.listdir("."):
        if ".log" in f or ".pckl" in f or ".traj" in f or "gfnff_topo" in f:
            os.remove(f)


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_static_Job():

    atoms = molecule("H2O")
    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["spin_multiplicity"] == 1
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["name"] == "xTB-Static"
    assert output["results"]["energy"] == pytest.approx(-137.96777587302995)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    job = StaticJob(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.96750578831137)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    job = StaticJob(method="GFN-FF").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN-FF"
    assert output["results"]["energy"] == pytest.approx(-8.912667188932252)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    atoms = bulk("Cu")
    job = StaticJob(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["results"]["energy"] == pytest.approx(-119.77643232313169)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_relax_Job():
    atoms = molecule("H2O")
    job = RelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["spin_multiplicity"] == 1
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["name"] == "xTB-Relax"
    assert output["results"]["energy"] == pytest.approx(-137.9764670127011)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    job = RelaxJob(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.9763496338962)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    job = RelaxJob(method="GFN-FF").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN-FF"
    assert output["results"]["energy"] == pytest.approx(-8.915974748299963)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    atoms = molecule("H2O")
    atoms.center(vacuum=5)
    atoms.pbc = True
    atoms[0].position += 0.01
    job = RelaxJob(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["results"]["energy"] == pytest.approx(-156.97441169886613)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_thermo_job():
    atoms = molecule("H2O")
    job = ThermoJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
