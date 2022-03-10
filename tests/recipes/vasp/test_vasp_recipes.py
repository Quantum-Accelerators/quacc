import os

from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

from quacc.recipes.vasp.core import DoubleRelaxMaker, RelaxMaker, StaticMaker
from quacc.recipes.vasp.qmof import QMOFMaker
from quacc.recipes.vasp.slabs import (
    BulkToSlabMaker,
    SlabRelaxMaker,
    SlabStaticMaker,
    SlabToAdsSlabMaker,
)


def teardown_module():
    for f in ["prerelax.log", "prerelax.traj"]:
        if os.path.exists(f):
            os.remove(f)


def test_static_maker():

    atoms = bulk("Cu") * (2, 2, 2)

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 2
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] == True
    assert output["name"] == "VASP-Static"

    job = StaticMaker(
        preset="BulkSet", name="test", swaps={"ncore": 2, "kpar": 4}
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4
    assert output["name"] == "test"

    job = StaticMaker(
        preset="QMOFSet", swaps={"ismear": 0, "sigma": 0.01, "nedos": None}
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["ismear"] == 0
    assert output["parameters"]["sigma"] == 0.01


def test_relax_maker():

    atoms = bulk("Cu") * (2, 2, 2)

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == False
    assert output["name"] == "VASP-Relax"

    job = RelaxMaker(preset="BulkSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"

    job = RelaxMaker(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["isif"] == 2


def test_doublerelax_maker():

    atoms = bulk("Cu") * (2, 2, 2)

    job = DoubleRelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == True
    assert output["name"] == "VASP-DoubleRelax"

    job = DoubleRelaxMaker(preset="BulkSet", name="test", swaps2={"nelmin": 6}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"

    job = DoubleRelaxMaker(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["isif"] == 2

    job = DoubleRelaxMaker(swaps1={"kpts": [1, 1, 1]}).make(atoms)
    responses = run_locally(job, ensure_success=True)


def test_slab_static_maker():
    atoms = bulk("Cu") * (2, 2, 2)

    job = SlabStaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] == True
    assert output["name"] == "VASP-SlabStatic"

    job = SlabStaticMaker(preset="SlabSet", name="test", swaps={"nelmin": 6}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"


def test_slab_relax_maker():
    atoms = bulk("Cu") * (2, 2, 2)

    job = SlabRelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] == False
    assert output["name"] == "VASP-SlabRelax"

    job = SlabRelaxMaker(preset="SlabSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["nelmin"] == 6
    assert output["name"] == "test"


def test_slab_flows():
    atoms = bulk("Cu") * (2, 2, 2)

    ### --------- Test BulkToSlabMaker --------- ###
    flow = BulkToSlabMaker().make(atoms)
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
    flow = BulkToSlabMaker(
        preset="SlabSet",
        name="test",
        slab_relax_maker=SlabRelaxMaker(swaps={"nelmin": 6}),
        slab_static_maker=SlabStaticMaker(swaps={"nelmin": 6}),
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

    ### --------- Test SlabToAdsSlabMaker --------- ###
    atoms = output2["atoms"]
    adsorbate = molecule("H2")

    flow = SlabToAdsSlabMaker().make(atoms, adsorbate)
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
    flow = SlabToAdsSlabMaker(preset="SlabSet", name="test", swaps={"nelmin": 6}).make(
        atoms, adsorbate
    )
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


def test_qmof_maker():
    atoms = bulk("Cu")
    job = QMOFMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["sigma"] == 0.01
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["name"] == "QMOF-Relax"

    job = QMOFMaker(preset="BulkSet", name="test", swaps={"nelmin": 6}).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["nelmin"] == 6
    assert output["parameters"]["sigma"] == 0.05
    assert output["name"] == "test"

    job = QMOFMaker(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    assert output["parameters"]["isif"] == 2

    atoms = bulk("Cu") * (8, 8, 8)
    job = QMOFMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
