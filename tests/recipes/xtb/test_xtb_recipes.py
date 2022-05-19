import os
from shutil import rmtree

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
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


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
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert len(output["results"]["frequencies"]) == 9
    assert len(output["results"]["true_frequencies"]) == 3
    assert output["results"]["true_frequencies"][-1] == pytest.approx(3526.945468014458)
    assert output["results"]["geometry"] == "nonlinear"
    assert output["results"]["energy"] == 0.0
    assert output["results"]["enthalpy"] == pytest.approx(0.637581401404518)
    assert output["results"]["entropy"] == pytest.approx(0.003942713004759747)
    assert output["results"]["gibbs_energy"] == pytest.approx(-0.5379384809646004)

    atoms = molecule("O2")
    job = ThermoJob(temperature=200, pressure=2.0).make(atoms, energy=-100.0)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert len(output["results"]["true_frequencies"]) == 1
    assert output["results"]["true_frequencies"][-1] == pytest.approx(1449.82397338371)
    assert output["results"]["geometry"] == "linear"
    assert output["results"]["pointgroup"] == "D*h"
    assert output["results"]["energy"] == -100.0
    assert output["results"]["enthalpy"] == pytest.approx(-99.84979574721257)
    assert output["results"]["entropy"] == pytest.approx(0.0038997333854535934)
    assert output["results"]["gibbs_energy"] == pytest.approx(-100.62974242430329)

    atoms = molecule("H")
    job = ThermoJob().make(atoms, energy=-1.0)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert output["results"]["geometry"] == "monatomic"
    assert len(output["results"]["true_frequencies"]) == 0
    assert output["results"]["energy"] == -1.0
    assert output["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
