from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.constraints import FixAtoms
from ase.md.npt import NPT
from ase.units import fs

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.md import md_job
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
    assert output["results"]["energy"] == pytest.approx(-0.04543069081693929)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01
    assert len(output["trajectory"]) == 30
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
    assert output["results"]["energy"] == pytest.approx(-0.004528885890177747)
    assert 0.01 < np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.03

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    c = FixAtoms(indices=[0, 1])
    atoms.set_constraint(c)
    output = relax_job(atoms, opt_params={"fmax": 0.03}, asap_cutoff=True)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(0.04996032884581858)


def test_md_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    old_positions = atoms.positions.copy()
    output = md_job(atoms)

    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 501
    assert output["name"] == "EMT Microcanonical"
    assert output["parameters_md"]["timestep"] == pytest.approx(1.0)
    assert output["trajectory_log"]["temperature"][-1] == pytest.approx(1575.886)
    assert output["trajectory_log"]["temperature"][0] == pytest.approx(0.0)
    assert output["trajectory_log"]["temperature"][1] == pytest.approx(759.680)
    assert output["trajectory_log"]["time"][10] == pytest.approx(0.01)
    assert atoms.positions == pytest.approx(old_positions)

    atoms = molecule("H2O")
    old_positions = atoms.positions.copy()

    rng = np.random.default_rng(seed=42)

    output = md_job(
        atoms,
        maxwell_boltzmann_params={"temperature": 1000, "rng": rng},
        md_params={"timestep": 0.5, "steps": 20},
    )

    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 21
    assert output["name"] == "EMT Microcanonical"
    assert output["parameters_md"]["timestep"] == pytest.approx(0.5)
    assert output["trajectory_log"]["temperature"][-1] == pytest.approx(1023.384)
    assert output["trajectory_log"]["temperature"][0] == pytest.approx(915.678)
    assert output["trajectory_log"]["temperature"][1] == pytest.approx(1060.650)
    assert output["trajectory_log"]["time"][10] == pytest.approx(0.005)
    assert atoms.positions == pytest.approx(old_positions)

    with pytest.raises(ValueError, match="Quacc does not support"):
        output = md_job(atoms, md_params={"dynamics_kwargs": {"trajectory": "md.traj"}})

    atoms = molecule("H2O", vacuum=10.0)
    old_positions = atoms.positions.copy()
    output = md_job(
        atoms,
        md_params={
            "timestep": 1.0,
            "dynamics": NPT,
            "dynamics_kwargs": {"temperature": 1000, "ttime": 50 * fs},
        },
    )

    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 500
    assert output["name"] == "EMT Microcanonical"


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
