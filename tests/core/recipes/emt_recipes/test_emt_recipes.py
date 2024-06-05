from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.constraints import FixAtoms
from ase.optimize import FIRE

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.info = {"test": "hello"}

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(0.07001766638245854)
    assert output["atoms"].info["test"] == "hello"
    assert output["atoms"].info.get("_id")

    output = static_job(atoms, asap_cutoff=True)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(0.11074520235398744)


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(-0.045446842063617154)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01
    assert len(output["trajectory"]) > 1
    assert output["atoms"] != output["input_atoms"]["atoms"]
    assert output["trajectory"][0] == output["input_atoms"]["atoms"]
    assert output["trajectory"][-1] == output["atoms"]
    assert (
        output["trajectory_results"][0]["energy"]
        > output["trajectory_results"][-1]["energy"]
    )
    assert output["trajectory_results"][-1]["energy"] == output["results"]["energy"]

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = relax_job(atoms, relax_cell=True)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    atoms = molecule("N2")
    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    output = relax_job(atoms, opt_params={"fmax": 0.03}, asap_cutoff=True)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(-0.004774645162642699)
    assert 0.01 < np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.03

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    c = FixAtoms(indices=[0, 1])
    atoms.set_constraint(c)
    output = relax_job(
        atoms, opt_params={"fmax": 0.03, "optimizer": FIRE}, asap_cutoff=True
    )
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(0.04996032884581858)


def test_slab_dynamic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    outputs = bulk_to_slabs_flow(atoms, run_static=False)
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 80
    assert outputs[1]["nsites"] == 96
    assert outputs[2]["nsites"] == 80
    assert outputs[3]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] is False for output in outputs]
    assert [output["name"] == "EMT Relax" for output in outputs]

    outputs = bulk_to_slabs_flow(
        atoms,
        run_static=False,
        job_params={"relax_job": {"opt_params": {"fmax": 1.0}, "asap_cutoff": True}},
    )
    assert len(outputs) == 4
    assert outputs[0]["nsites"] == 80
    assert outputs[1]["nsites"] == 96
    assert outputs[2]["nsites"] == 80
    assert outputs[3]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] is True for output in outputs]


def test_customizer():
    atoms = bulk("Cu")
    results = bulk_to_slabs_flow(
        atoms, job_params={"static_job": {"asap_cutoff": True}}
    )
    for result in results:
        assert result["parameters"]["asap_cutoff"] is True


def test_customizer_v2():
    atoms = bulk("Cu")
    results = bulk_to_slabs_flow(atoms, job_params={"relax_job": {"asap_cutoff": True}})
    for result in results:
        assert result["parameters"]["asap_cutoff"] is False
        assert Path(result["dir_name"], "quacc_results.json.gz").exists()


def test_all_customizers():
    atoms = bulk("Cu")
    results = bulk_to_slabs_flow(atoms, job_params={"all": {"asap_cutoff": True}})
    for result in results:
        assert result["parameters"]["asap_cutoff"] is True


def test_all_customizers_v2():
    atoms = bulk("Cu")
    results = bulk_to_slabs_flow(
        atoms,
        job_params={"all": {"asap_cutoff": True}, "static_job": {"asap_cutoff": False}},
    )
    for result in results:
        assert result["parameters"]["asap_cutoff"] is False
