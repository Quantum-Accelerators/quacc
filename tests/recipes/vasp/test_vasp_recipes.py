import pytest

from ase.build import bulk, molecule
from jobflow.managers.local import run_locally

from quacc.schemas.calc import summarize_run as calc_summarize_run
from quacc.recipes.vasp.core import RelaxMaker, StaticMaker
from quacc.recipes.vasp.slabs import (
    BulkToSlabMaker,
    SlabRelaxMaker,
    SlabStaticMaker,
    SlabToAdsSlabMaker,
)
from quacc.util.json import jsonify, unjsonify


def mock_summarize_run(atoms, **kwargs):
    # Instead of running the VASP-specific summarize_run(), we mock it with the
    # general calculator schema which does not require VASP files to be
    # in the working directory and will work with pytest.

    prep_next_run = kwargs.get("prep_next_run", True)
    additioanl_fields = kwargs.get("additional_fields", None)
    output = calc_summarize_run(
        atoms, prep_next_run=prep_next_run, additional_fields=additioanl_fields
    )
    return output


@pytest.fixture(autouse=True)
def patch_summarize_run(monkeypatch):
    # Monkeypatch the summarize_run() function so that we aren't relying on real
    # VASP files to be in the working directory during the test. Note that even though
    # summarize_run() is a function in the quacc.schemas.vasp module, we modify it
    # only in quacc.recipes.vasp.core/.slabs because otherwise it will not work properly.
    monkeypatch.setattr("quacc.recipes.vasp.core.summarize_run", mock_summarize_run)
    monkeypatch.setattr("quacc.recipes.vasp.slabs.summarize_run", mock_summarize_run)


def test_static_maker():

    atoms = bulk("Cu") * (2, 2, 2)
    atoms_json = jsonify(atoms)

    job = StaticMaker().make(atoms_json)
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 2
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] == True
    assert output["name"] == "Static"

    job = StaticMaker(preset="BulkRelaxSet", ncore=2, kpar=4, name="test").make(
        atoms_json
    )
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4
    assert output["name"] == "test"


def test_relax_maker():

    atoms = bulk("Cu") * (2, 2, 2)
    atoms_json = jsonify(atoms)

    job = RelaxMaker().make(atoms_json)
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isif"] == 3
    assert output["parameters"]["lwave"] == False
    assert output["name"] == "Relax"

    job = RelaxMaker(preset="BulkRelaxSet", ncore=2, kpar=4, name="test").make(
        atoms_json
    )
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 650
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4
    assert output["name"] == "test"


def test_slab_static_maker():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms_json = jsonify(atoms)

    job = SlabStaticMaker().make(atoms_json)
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["idipol"] == 3
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lvhar"] == True
    assert output["name"] == "SlabStatic"

    job = SlabStaticMaker(preset="SlabRelaxSet", ncore=2, kpar=4, name="test").make(
        atoms_json
    )
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4
    assert output["name"] == "test"


def test_slab_relax_maker():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms_json = jsonify(atoms)

    job = SlabRelaxMaker().make(atoms_json)
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["isif"] == 2
    assert output["parameters"]["nsw"] > 0
    assert output["parameters"]["isym"] == 0
    assert output["parameters"]["lwave"] == False
    assert output["name"] == "SlabRelax"

    job = SlabRelaxMaker(preset="SlabRelaxSet", ncore=2, kpar=4, name="test").make(
        atoms_json
    )
    responses = run_locally(job)
    output = responses[job.uuid][1].output
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["ncore"] == 2
    assert output["parameters"]["kpar"] == 4
    assert output["name"] == "test"


def test_slab_flows():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms_json = jsonify(atoms)

    ## Test BulkToSlabMaker
    flow = BulkToSlabMaker().make(atoms_json)
    responses = run_locally(flow)

    assert len(responses) == 9
    uuids = list(responses.keys())

    # First job is a dummy job to make slabs and should have no output
    output0 = responses[uuids[0]][1].output
    assert output0 is None

    # Subsequent jobs should be alternating relaxations and statics
    output1 = responses[uuids[1]][1].output
    assert output1["nsites"] > len(atoms)
    assert output1["parameters"]["isif"] == 2
    assert output1["name"] == "BulkToSlab_SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["nsites"] == output1["nsites"]
    assert output2["parameters"]["nsw"] == 0
    assert output2["name"] == "BulkToSlab_SlabStatic"

    # Now try with kwargs
    flow = BulkToSlabMaker(preset="SlabRelaxSet", ncore=2, kpar=4, name="test").make(
        atoms_json
    )
    responses = run_locally(flow)

    assert len(responses) == 9
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    assert output0 is None

    output1 = responses[uuids[1]][1].output
    assert output1["parameters"]["isif"] == 2
    assert output1["parameters"]["ncore"] == 2
    assert output1["parameters"]["kpar"] == 4
    assert output1["parameters"]["encut"] == 450
    assert output1["name"] == "test_SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["parameters"]["nsw"] == 0
    assert output2["parameters"]["ncore"] == 2
    assert output2["parameters"]["kpar"] == 4
    assert output2["parameters"]["encut"] == 450
    assert output2["name"] == "test_SlabStatic"

    ## Test SlabtoAdsSlabMaker using one of the slabs from above
    atoms_json = output2["atoms"]
    adsorbate = molecule("H2")
    adsorbate_json = jsonify(adsorbate)

    flow = SlabToAdsSlabMaker().make(atoms_json, adsorbate_json)
    responses = run_locally(flow)

    assert len(responses) == 11
    uuids = list(responses.keys())

    # First job is a dummy job to make slabs and should have no output
    output0 = responses[uuids[0]][1].output
    assert output0 is None

    # Subsequent jobs should be alternating relaxations and statics
    output1 = responses[uuids[1]][1].output
    assert output1["nsites"] == len(unjsonify(output2["atoms"])) + 2
    assert output1["parameters"]["isif"] == 2
    assert output1["name"] == "SlabToAdsSlab_SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["nsites"] == output1["nsites"]
    assert output2["parameters"]["nsw"] == 0
    assert output2["name"] == "SlabToAdsSlab_SlabStatic"

    # Now try with kwargs
    flow = SlabToAdsSlabMaker(preset="SlabRelaxSet", ncore=2, kpar=4, name="test").make(
        atoms_json, adsorbate_json
    )
    responses = run_locally(flow)

    assert len(responses) == 11
    uuids = list(responses.keys())

    output0 = responses[uuids[0]][1].output
    assert output0 is None

    output1 = responses[uuids[1]][1].output
    assert output1["parameters"]["isif"] == 2
    assert output1["parameters"]["ncore"] == 2
    assert output1["parameters"]["kpar"] == 4
    assert output1["parameters"]["encut"] == 450
    assert output1["name"] == "test_SlabRelax"

    output2 = responses[uuids[2]][1].output
    assert output2["parameters"]["nsw"] == 0
    assert output2["parameters"]["ncore"] == 2
    assert output2["parameters"]["kpar"] == 4
    assert output2["parameters"]["encut"] == 450
    assert output2["name"] == "test_SlabStatic"
