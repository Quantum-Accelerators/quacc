from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

from quacc.recipes.vasp.core import RelaxMaker, StaticMaker
from quacc.recipes.vasp.slabs import (
    BulkToSlabMaker,
    SlabRelaxMaker,
    SlabStaticMaker,
    SlabToAdsSlabMaker,
)


def test_static_maker():

    atoms = bulk("Cu") * (2, 2, 2)

    job = StaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["isym"] != 2:
        raise AssertionError
    if output["parameters"]["nsw"] != 0:
        raise AssertionError
    if output["parameters"]["lwave"] != True:
        raise AssertionError
    if output["name"] != "VASP-Static":
        raise AssertionError

    job = StaticMaker(
        preset="BulkRelaxSet", name="test", swaps={"ncore": 2, "kpar": 4}
    ).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["parameters"]["encut"] != 650:
        raise AssertionError
    if output["parameters"]["ncore"] != 2:
        raise AssertionError
    if output["parameters"]["kpar"] != 4:
        raise AssertionError
    if output["name"] != "test":
        raise AssertionError


def test_relax_maker():

    atoms = bulk("Cu") * (2, 2, 2)

    job = RelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["isym"] != 0:
        raise AssertionError
    if output["parameters"]["nsw"] <= 0:
        raise AssertionError
    if output["parameters"]["isif"] != 3:
        raise AssertionError
    if output["parameters"]["lwave"] != False:
        raise AssertionError
    if output["name"] != "VASP-Relax":
        raise AssertionError

    job = RelaxMaker(preset="BulkRelaxSet", name="test", swaps={"nelmin": 6}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["parameters"]["encut"] != 650:
        raise AssertionError
    if output["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output["name"] != "test":
        raise AssertionError

    job = RelaxMaker(volume_relax=False).make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["parameters"]["isif"] != 2:
        raise AssertionError


def test_slab_static_maker():
    atoms = bulk("Cu") * (2, 2, 2)

    job = SlabStaticMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["idipol"] != 3:
        raise AssertionError
    if output["parameters"]["nsw"] != 0:
        raise AssertionError
    if output["parameters"]["lvhar"] != True:
        raise AssertionError
    if output["name"] != "VASP-SlabStatic":
        raise AssertionError

    job = SlabStaticMaker(preset="SlabRelaxSet", name="test", swaps={"nelmin": 6}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["parameters"]["encut"] != 450:
        raise AssertionError
    if output["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output["name"] != "test":
        raise AssertionError


def test_slab_relax_maker():
    atoms = bulk("Cu") * (2, 2, 2)

    job = SlabRelaxMaker().make(atoms)
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["nsites"] != len(atoms):
        raise AssertionError
    if output["parameters"]["isif"] != 2:
        raise AssertionError
    if output["parameters"]["nsw"] <= 0:
        raise AssertionError
    if output["parameters"]["isym"] != 0:
        raise AssertionError
    if output["parameters"]["lwave"] != False:
        raise AssertionError
    if output["name"] != "VASP-SlabRelax":
        raise AssertionError

    job = SlabRelaxMaker(preset="SlabRelaxSet", name="test", swaps={"nelmin": 6}).make(
        atoms
    )
    responses = run_locally(job, ensure_success=True)
    output = responses[job.uuid][1].output
    if output["parameters"]["encut"] != 450:
        raise AssertionError
    if output["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output["name"] != "test":
        raise AssertionError


def test_slab_flows():
    atoms = bulk("Cu") * (2, 2, 2)

    ### --------- Test BulkToSlabMaker --------- ###
    flow = BulkToSlabMaker().make(atoms)
    responses = run_locally(flow, ensure_success=True)

    if len(responses) != 9:
        raise AssertionError
    uuids = list(responses.keys())

    # First job is a dummy job to make slabs and should have no output
    output0 = responses[uuids[0]][1].output
    if output0 is not None:
        raise AssertionError

    output1 = responses[uuids[1]][1].output
    if output1["nsites"] <= len(atoms):
        raise AssertionError
    if output1["parameters"]["isif"] != 2:
        raise AssertionError
    if output1["name"] != "VASP-SlabRelax":
        raise AssertionError

    output2 = responses[uuids[2]][1].output
    if output2["nsites"] != output1["nsites"]:
        raise AssertionError
    if output2["parameters"]["nsw"] != 0:
        raise AssertionError
    if output2["name"] != "VASP-SlabStatic":
        raise AssertionError

    # Now try with kwargs
    flow = BulkToSlabMaker(
        preset="SlabRelaxSet",
        name="test",
        slab_relax_maker=SlabRelaxMaker(swaps={"nelmin": 6}),
        slab_static_maker=SlabStaticMaker(swaps={"nelmin": 6}),
    ).make(atoms)
    responses = run_locally(flow, ensure_success=True)

    if len(responses) != 9:
        raise AssertionError
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    if output0 is not None:
        raise AssertionError

    output1 = responses[uuids[1]][1].output
    if output1["parameters"]["isif"] != 2:
        raise AssertionError
    if output1["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output1["parameters"]["encut"] != 450:
        raise AssertionError
    if output1["name"] != "VASP-SlabRelax":
        raise AssertionError

    output2 = responses[uuids[2]][1].output
    if output2["parameters"]["nsw"] != 0:
        raise AssertionError
    if output2["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output2["parameters"]["encut"] != 450:
        raise AssertionError
    if output2["name"] != "VASP-SlabStatic":
        raise AssertionError

    ### --------- Test SlabToAdsSlabMaker --------- ###
    atoms = output2["atoms"]
    adsorbate = molecule("H2")

    flow = SlabToAdsSlabMaker().make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)

    if len(responses) != 11:
        raise AssertionError
    uuids = list(responses.keys())

    # First job is a dummy job to make slabs and should have no output
    output0 = responses[uuids[0]][1].output
    if output0 is not None:
        raise AssertionError

    # Subsequent jobs should be alternating relaxations and statics
    output1 = responses[uuids[1]][1].output
    if output1["nsites"] != len(output2["atoms"]) + 2:
        raise AssertionError
    if output1["parameters"]["isif"] != 2:
        raise AssertionError
    if output1["name"] != "VASP-SlabRelax":
        raise AssertionError

    output2 = responses[uuids[2]][1].output
    if output2["nsites"] != output1["nsites"]:
        raise AssertionError
    if output2["parameters"]["nsw"] != 0:
        raise AssertionError
    if output2["name"] != "VASP-SlabStatic":
        raise AssertionError

    # Now try with kwargs
    flow = SlabToAdsSlabMaker(
        preset="SlabRelaxSet", name="test", swaps={"nelmin": 6}
    ).make(atoms, adsorbate)
    responses = run_locally(flow, ensure_success=True)

    if len(responses) != 11:
        raise AssertionError
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    if output0 is not None:
        raise AssertionError

    output1 = responses[uuids[1]][1].output
    if output1["parameters"]["isif"] != 2:
        raise AssertionError
    if output1["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output1["parameters"]["encut"] != 450:
        raise AssertionError
    if output1["name"] != "VASP-SlabRelax":
        raise AssertionError

    output2 = responses[uuids[2]][1].output
    if output2["parameters"]["nsw"] != 0:
        raise AssertionError
    if output2["parameters"]["nelmin"] != 6:
        raise AssertionError
    if output2["parameters"]["encut"] != 450:
        raise AssertionError
    if output2["name"] != "VASP-SlabStatic":
        raise AssertionError
