import os

import pytest
from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

from quacc.recipes.vasp.core import DoubleRelaxJob, RelaxJob, StaticJob
from quacc.recipes.vasp.qmof import QMOFRelaxJob
from quacc.recipes.vasp.slabs import (
    BulkToAdsorbatesFlow,
    BulkToSlabsJob,
    SlabRelaxJob,
    SlabStaticJob,
    SlabToAdsorbatesJob,
)


def teardown_module():
    for f in ["prerelax.log", "prerelax.traj"]:
        if os.path.exists(f):
            os.remove(f)


def test_static_job():

    atoms = bulk("Cu") * (2, 2, 2)

    job = StaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 2
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] == True
    assert output["name"] == "VASP-Static"

    job = StaticJob(preset="BulkSet", name="test", swaps={"ncore": 2, "kpar": 4}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4
    assert output["name"] == "test"

    job = StaticJob(
        preset="QMOFSet", swaps={"ismear": 0, "sigma": 0.01, "nedos": None}
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.01


def test_relax_job():

    atoms = bulk("Cu") * (2, 2, 2)

    job = RelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == False
    assert output["name"] == "VASP-Relax"

    job = RelaxJob(preset="BulkSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"

    job = RelaxJob(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["isif"] == 2


def test_doublerelax_job():

    atoms = bulk("Cu") * (2, 2, 2)

    job = DoubleRelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == True
    assert output["name"] == "VASP-DoubleRelax"

    job = DoubleRelaxJob(preset="BulkSet", name="test", swaps2={"nelmin": 6}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"

    job = DoubleRelaxJob(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["isif"] == 2

    job = DoubleRelaxJob(swaps1={"kpts": [1, 1, 1]}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output


def test_slab_static_job():
    atoms = bulk("Cu") * (2, 2, 2)

    job = SlabStaticJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] == True
    assert output["name"] == "VASP-SlabStatic"

    job = SlabStaticJob(preset="SlabSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"


def test_slab_relax_job():
    atoms = bulk("Cu") * (2, 2, 2)

    job = SlabRelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] == False
    assert output["name"] == "VASP-SlabRelax"

    job = SlabRelaxJob(preset="SlabSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"


def test_slab_dynamic_jobs():
    atoms = bulk("Cu") * (2, 2, 2)

    ### --------- Test BulkToSlabsJob --------- ###
    flow = BulkToSlabsJob().make(atoms)
    responses = run_locally(flow, ensure_success=True)

    assert len(responses) == 9
    uuids = list(responses.keys())

    # First job is a dummy job to make slabs and should have no output
    output0 = responses[uuids[0]][1].output
    assert output0 is None

    output1 = responses[uuids[1]][1].output
    assert output1["nsites"] > len(atoms)
    assert output1["parameters"]["isif"] == 2
    assert output1["name"] == "VASP-SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["nsites"] == output1["nsites"]
    assert output2["parameters"]["nsw"] == 0
    assert output2["name"] == "VASP-SlabStatic"

    # Now try with kwargs
    flow = BulkToSlabsJob(
        name="test",
        slab_relax_job=SlabRelaxJob(preset="SlabSet", swaps={"nelmin": 6}),
        slab_static_job=SlabStaticJob(preset="SlabSet", swaps={"nelmin": 6}),
    ).make(atoms)
    responses = run_locally(flow, ensure_success=True)

    assert len(responses) == 9
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    assert output0 is None

    output1 = responses[uuids[1]][1].output
    assert output1["parameters"]["isif"] == 2
    assert output1["parameters"]["nelmin"] == 6
    assert output1["parameters"]["encut"] == 450
    assert output1["name"] == "VASP-SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["parameters"]["nsw"] == 0
    assert output2["parameters"]["nelmin"] == 6
    assert output2["parameters"]["encut"] == 450
    assert output2["name"] == "VASP-SlabStatic"

    ### --------- Test SlabToAdsorbatesJob --------- ###
    atoms = output2["atoms"]
    adsorbate = molecule("H2")

    flow = SlabToAdsorbatesJob().make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)

    assert len(responses) == 11
    uuids = list(responses.keys())

    # First job is a dummy job to make slabs and should have no output
    output0 = responses[uuids[0]][1].output
    assert output0 is None

    # Subsequent jobs should be alternating relaxations and statics
    output1 = responses[uuids[1]][1].output
    assert output1["nsites"] == len(output2["atoms"]) + 2
    assert output1["parameters"]["isif"] == 2
    assert output1["name"] == "VASP-SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["nsites"] == output1["nsites"]
    assert output2["parameters"]["nsw"] == 0
    assert output2["name"] == "VASP-SlabStatic"

    # Now try with kwargs
    flow = SlabToAdsorbatesJob(
        name="test",
        slab_relax_job=SlabRelaxJob(preset="SlabSet", swaps={"nelmin": 6}),
        slab_static_job=SlabStaticJob(preset="SlabSet", swaps={"nelmin": 6}),
    ).make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)

    assert len(responses) == 11
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    assert output0 is None

    output1 = responses[uuids[1]][1].output
    assert output1["parameters"]["isif"] == 2
    assert output1["parameters"]["nelmin"] == 6
    assert output1["parameters"]["encut"] == 450
    assert output1["name"] == "VASP-SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["parameters"]["nsw"] == 0
    assert output2["parameters"]["nelmin"] == 6
    assert output2["parameters"]["encut"] == 450
    assert output2["name"] == "VASP-SlabStatic"


def test_slab_flows():
    # TODO: This could use some more detailed tests

    atoms = bulk("Cu")
    adsorbate = molecule("H2O")
    flow = BulkToAdsorbatesFlow().make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)
    uuids = list(responses.keys())
    assert len(responses) == 27

    output0 = responses[uuids[0]][1].output
    assert output0["parameters"]["ediffg"] == -0.02

    output1 = responses[uuids[1]][1].output
    assert output1["parameters"]["nsw"] == 0

    output_final = responses[uuids[-1]][1].output
    assert output_final["parameters"]["nsw"] == 0

    flow = BulkToAdsorbatesFlow(stable_slab=False).make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)
    assert len(responses) == 78

    flow = BulkToAdsorbatesFlow(
        bulk_relax_job=None, bulk_static_job=None, stable_slab=False
    ).make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)
    assert len(responses) == 76

    with pytest.raises(ValueError):
        flow = BulkToAdsorbatesFlow(bulk_relax_job=None, bulk_static_job=None).make(
            atoms, adsorbate
        )


def test_qmof():
    atoms = bulk("Cu")
    job = QMOFRelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["sigma"] == 0.01
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["name"] == "QMOF-Relax"

    job = QMOFRelaxJob(preset="BulkSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["nelmin"] == 6
    assert output["parameters"]["sigma"] == 0.05
    assert output["name"] == "test"

    job = QMOFRelaxJob(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["isif"] == 2

    atoms = bulk("Cu") * (8, 8, 8)
    job = QMOFRelaxJob().make(atoms)
    responses = run_locally(job, ensure_success=True)
