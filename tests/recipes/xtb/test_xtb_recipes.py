import os

import pytest
from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

try:
    from xtb.ase.calculator import XTB
except ModuleNotFoundError:
    XTB = None
from quacc.recipes.xtb.core import RelaxMaker, StaticMaker


def teardown_module():
    for f in os.listdir("."):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f == "gfnff_topo"
        ):
            os.remove(f)


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_static_maker():

    atoms = molecule("H2O")
    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["spin_multiplicity"] == 1
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["name"] == "xTB-Static"
    assert output["results"]["energy"] == pytest.approx(-137.96777587302995)

    job = StaticMaker(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.96750578831137)

    job = StaticMaker(method="GFN-FF").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN-FF"
    assert output["results"]["energy"] == pytest.approx(-8.912667188932252)

    atoms = bulk("Cu")
    job = StaticMaker(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["results"]["energy"] == pytest.approx(-119.77643232313169)


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_relax_maker():
    atoms = molecule("H2O")
    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["spin_multiplicity"] == 1
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["name"] == "xTB-Relax"
    assert output["results"]["energy"] == pytest.approx(-137.9764670127011)

    job = RelaxMaker(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.9763496338962)

    job = RelaxMaker(method="GFN-FF").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["method"] == "GFN-FF"
    assert output["results"]["energy"] == pytest.approx(-8.915963671549369)

    atoms = bulk("Cu")
    job = RelaxMaker(method="GFN1-xTB").make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["results"]["energy"] == pytest.approx(-119.77643232313169)
