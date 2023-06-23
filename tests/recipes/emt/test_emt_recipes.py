import os
from shutil import rmtree

import jobflow as jf
import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.constraints import FixAtoms
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.jobflow.slabs import BulkToSlabsFlow as JFBulkToSlabsFlow
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


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
        if "quacc-tmp-" in f or f == "tmp_dir":
            rmtree(f)


def test_static_Job():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(0.07001766638245854)

    output = static_job(atoms, calc_kwargs={"asap_cutoff": True})
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(0.11074520235398744)


def test_relax_Job():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = relax_job(atoms, relax_cell=False)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(-0.04543069081693929)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    atoms = molecule("N2")
    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    output = relax_job(
        atoms,
        relax_cell=False,
        opt_swaps={"fmax": 0.03},
        calc_kwargs={"asap_cutoff": True},
    )
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(-0.004528885890177747)
    assert 0.01 < np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.03

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    c = FixAtoms(indices=[0, 1])
    atoms.set_constraint(c)
    output = relax_job(
        atoms,
        relax_cell=False,
        opt_swaps={"fmax": 0.03},
        calc_kwargs={"asap_cutoff": True},
    )
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(0.04996032884581858)


def test_slab_dynamic_jobs():
    atoms = bulk("Cu")

    outputs = bulk_to_slabs_flow(atoms, slab_static_electron=None)
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 80
    assert outputs[1]["nsites"] == 96
    assert outputs[2]["nsites"] == 80
    assert outputs[3]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] is False for output in outputs]
    assert [output["name"] == "EMT Relax" for output in outputs]

    outputs = bulk_to_slabs_flow(
        atoms,
        slab_static_electron=None,
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_kwargs": {"asap_cutoff": True},
            "relax_cell": False,
        },
    )
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 80
    assert outputs[1]["nsites"] == 96
    assert outputs[2]["nsites"] == 80
    assert outputs[3]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] is True for output in outputs]

    outputs = bulk_to_slabs_flow(
        atoms,
        slabgen_kwargs={"max_slabs": 2},
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_kwargs": {"asap_cutoff": True},
            "relax_cell": False,
        },
    )
    assert len(outputs) == 2
    assert outputs[0]["nsites"] == 64
    assert outputs[1]["nsites"] == 80
    assert [output["parameters"]["asap_cutoff"] is False for output in outputs]


def test_jf_slab_dynamic_jobs():
    store = jf.JobStore(MemoryStore())

    atoms = bulk("Cu")

    flow = JFBulkToSlabsFlow(slab_static_job=None).make(atoms)
    jf.run_locally(flow, store=store, create_folders=True, ensure_success=True)

    flow = JFBulkToSlabsFlow(
        slab_static_job=None,
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_kwargs": {"asap_cutoff": True},
            "relax_cell": False,
        },
    ).make(atoms)
    jf.run_locally(flow, store=store, create_folders=True, ensure_success=True)

    flow = JFBulkToSlabsFlow(
        slab_relax_kwargs={
            "opt_swaps": {"fmax": 1.0},
            "calc_kwargs": {"asap_cutoff": True},
            "relax_cell": False,
        },
    ).make(atoms, slabgen_kwargs={"max_slabs": 2})
    responses = jf.run_locally(
        flow, store=store, create_folders=True, ensure_success=True
    )

    assert len(responses) == 5
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    assert "generated_slabs" in output0
    assert len(output0["generated_slabs"][0]) == 64

    output1 = responses[uuids[1]][1].output
    assert output1["nsites"] == 64
    assert output1["parameters"]["asap_cutoff"] is True
    assert output1["name"] == "EMT Relax"

    output2 = responses[uuids[-1]][1].output
    assert output2["nsites"] == 80
    assert output2["parameters"]["asap_cutoff"] is False
    assert output2["name"] == "EMT Static"
