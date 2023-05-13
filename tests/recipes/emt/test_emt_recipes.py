import os
from shutil import rmtree

import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job, static_job


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or ".gz" in f
        ):
            os.remove(f)
        if "quacc-tmp" in f or f == "tmp_dir":
            rmtree(f)


def test_static_Job():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] == False
    assert output["results"]["energy"] == pytest.approx(0.07001766638245854)

    output = static_job(atoms, emt_kwargs={"asap_cutoff": True})
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] == True
    assert output["results"]["energy"] == pytest.approx(0.11074520235398744)


def test_relax_Job():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] == False
    assert output["results"]["energy"] == pytest.approx(-0.04543069081693929)

    output = relax_job(atoms, fmax=0.03, emt_kwargs={"asap_cutoff": True})
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] == True
    assert output["results"]["energy"] == pytest.approx(-0.004528885890177747)
