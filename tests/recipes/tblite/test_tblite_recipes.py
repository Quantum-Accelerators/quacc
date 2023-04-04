import os

import numpy as np
import pytest
from ase.build import molecule
from jobflow.managers.local import run_locally

try:
    from tblite.ase import TBLite
except ModuleNotFoundError:
    TBLite = None
from quacc.recipes.tblite.core import RelaxJob, StaticJob


def teardown_module():
    for f in os.listdir("."):
        if ".log" in f or ".pckl" in f or ".traj" in f or "gfnff_topo" in f:
            os.remove(f)


@pytest.mark.skipif(
    TBLite is None,
    reason="tblite must be installed. Try conda install -c conda-forge tblite",
)
def test_static_Job():
    atoms = molecule("H2O")
    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["spin_multiplicity"] == 1
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["name"] == "tblite-Static"
    assert output["results"]["energy"] == pytest.approx(-137.9745971808714)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    job = StaticJob(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.96750578831137)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    TBLite is None,
    reason="tblite must be installed. Try conda install -c conda-forge tblite",
)
def test_relax_Job():
    atoms = molecule("H2O")
    job = RelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["spin_multiplicity"] == 1
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["name"] == "tblite-Relax"
    assert output["results"]["energy"] == pytest.approx(-137.9833607328507)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
