import os
from shutil import rmtree

import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow


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


def test_slab_dynamic_jobs():
    atoms = bulk("Cu")

    with pytest.raises(ValueError):
        BulkToSlabsFlow(relax_electron=None, static_electron=None).run(atoms)

    outputs = BulkToSlabsFlow(relax_electron=None).run(atoms)
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 80
    assert outputs[1]["nsites"] == 96
    assert outputs[2]["nsites"] == 80
    assert outputs[3]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] == False for output in outputs]

    outputs = BulkToSlabsFlow(
        static_electron=None,
        relax_kwargs={"fmax": 1.0, "emt_kwargs": {"asap_cutoff": True}},
    ).run(atoms)
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 80
    assert outputs[1]["nsites"] == 96
    assert outputs[2]["nsites"] == 80
    assert outputs[3]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] == True for output in outputs]

    outputs = BulkToSlabsFlow(
        relax_kwargs={"fmax": 1.0, "emt_kwargs": {"asap_cutoff": True}},
    ).run(atoms, slabgen_kwargs={"max_slabs": 2})
    assert len(outputs) == 2
    assert outputs[0]["nsites"] == 64
    assert outputs[1]["nsites"] == 80
    assert [output["parameters"]["asap_cutoff"] == False for output in outputs]
