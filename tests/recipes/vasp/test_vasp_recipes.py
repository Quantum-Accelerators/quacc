import pytest

from ase.build import bulk
from jobflow.managers.local import run_locally

from quacc.schemas.calc import summarize_run as calc_summarize_run
from quacc.recipes.vasp.core import RelaxMaker, StaticMaker
from quacc.recipes.vasp.slabs import (
    BulkToSlabMaker,
    SlabRelaxMaker,
    SlabStaticMaker,
    SlabToAdsSlabMaker,
)
from quacc.util.json import jsonify


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


def test_bulk_to_slab():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms_json = jsonify(atoms)

    job = BulkToSlabMaker().make(atoms_json)
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
