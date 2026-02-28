from __future__ import annotations

from logging import getLogger

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.constraints import FixAtoms
from ase.md.npt import NPT
from ase.optimize import FIRE
from ase.units import fs

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.md import md_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

LOGGER = getLogger(__name__)
LOGGER.propagate = True


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.info = {"test": "hello"}

    output = static_job(atoms)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert output["results"]["energy"] == pytest.approx(0.07001766638245854)
    assert output["atoms"].info["test"] == "hello"
    assert output["atoms"].info.get("_id")

    output = static_job(atoms, asap_cutoff=True)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(0.11074520235398744)


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = relax_job(atoms)
    assert output["structure_metadata"]["nsites"] == len(atoms)
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
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    atoms = molecule("N2")
    output = relax_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is False
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    output = relax_job(atoms, opt_params={"fmax": 0.03}, asap_cutoff=True)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert output["parameters"]["asap_cutoff"] is True
    assert output["results"]["energy"] == pytest.approx(-0.004774645162642699)
    assert 0.01 < np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.03

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    c = FixAtoms(indices=[0, 1])
    atoms.set_constraint(c)
    output_fire = relax_job(
        atoms, opt_params={"fmax": 0.03, "optimizer": FIRE}, asap_cutoff=True
    )
    assert output_fire["structure_metadata"]["nsites"] == len(atoms)
    assert output_fire["parameters"]["asap_cutoff"] is True
    assert output_fire["results"]["energy"] == pytest.approx(0.04996032884581858)

    # Add a test that passes through kwargs to the FrechetCellFilter
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]
    c = FixAtoms(indices=[0, 1])
    atoms.set_constraint(c)
    output_fire_pressure = relax_job(
        atoms,
        relax_cell=True,
        opt_params={
            "fmax": 0.03,
            "optimizer": FIRE,
            "filter_kwargs": {"scalar_pressure": 0.01},
        },
        asap_cutoff=True,
    )

    # Check that volume is smaller with a pressure applied than without!
    assert (
        output_fire_pressure["structure_metadata"]["volume"]
        < output_fire["structure_metadata"]["volume"]
    )
    assert output_fire_pressure["structure_metadata"]["nsites"] == len(atoms)
    assert output_fire_pressure["parameters"]["asap_cutoff"] is True


def test_md_job1():
    atoms = molecule("H2O")
    old_positions = atoms.positions.copy()
    output = md_job(atoms, steps=500)
    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 501
    assert output["name"] == "EMT MD"
    assert output["parameters_md"]["timestep"] == pytest.approx(1.0 * fs)
    assert output["trajectory_log"][-1]["temperature"] == pytest.approx(1575.886)
    assert output["trajectory_log"][0]["temperature"] == pytest.approx(0.0)
    assert output["trajectory_log"][1]["temperature"] == pytest.approx(759.680)
    assert output["trajectory_log"][10]["time"] == pytest.approx(10 * fs)
    assert atoms.positions == pytest.approx(old_positions)


def test_md_job2():
    atoms = molecule("H2O")
    old_positions = atoms.positions.copy()

    output = md_job(
        atoms,
        timestep_fs=0.5,
        steps=20,
        md_params={
            "maxwell_boltzmann_kwargs": {
                "temperature_K": 1000,
                "rng": np.random.default_rng(seed=42),
            },
            "set_com_stationary": True,
            "set_zero_rotation": True,
        },
    )
    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 21
    assert output["name"] == "EMT MD"
    assert output["parameters_md"]["timestep"] == pytest.approx(0.5 * fs)
    assert output["trajectory_log"][-1]["temperature"] == pytest.approx(1023.384)
    assert output["trajectory_log"][0]["temperature"] == pytest.approx(915.678)
    assert output["trajectory_log"][1]["temperature"] == pytest.approx(1060.650)
    assert output["trajectory_log"][10]["time"] == pytest.approx(10 * 0.5 * fs)
    assert atoms.positions == pytest.approx(old_positions)


def test_md_job3():
    atoms = molecule("H2O", vacuum=10.0)
    output = md_job(
        atoms,
        dynamics=NPT,
        timestep_fs=1.0,
        temperature_K=1000,
        steps=500,
        md_params={"dynamics_kwargs": {"ttime": 50 * fs}},
    )
    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 501
    assert output["name"] == "EMT MD"
    assert output["trajectory_log"][1]["temperature"] == pytest.approx(759.8829)
    assert output["trajectory_results"][-1]["energy"] == pytest.approx(2.0363759)


def test_md_job_error():
    atoms = molecule("H2O")
    with pytest.raises(ValueError, match="Quacc does not support"):
        md_job(atoms, md_params={"dynamics_kwargs": {"trajectory": "md.traj"}})


def test_slab_dynamic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    outputs = bulk_to_slabs_flow(atoms, run_static=False)
    assert len(outputs) == 4
    assert outputs[0]["structure_metadata"]["nsites"] == 80
    assert outputs[1]["structure_metadata"]["nsites"] == 96
    assert outputs[2]["structure_metadata"]["nsites"] == 80
    assert outputs[3]["structure_metadata"]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] is False for output in outputs]
    assert [output["name"] == "EMT Relax" for output in outputs]

    outputs = bulk_to_slabs_flow(
        atoms,
        run_static=False,
        job_params={"relax_job": {"opt_params": {"fmax": 1.0}, "asap_cutoff": True}},
    )
    assert len(outputs) == 4
    assert outputs[0]["structure_metadata"]["nsites"] == 80
    assert outputs[1]["structure_metadata"]["nsites"] == 96
    assert outputs[2]["structure_metadata"]["nsites"] == 80
    assert outputs[3]["structure_metadata"]["nsites"] == 64
    assert [output["parameters"]["asap_cutoff"] is True for output in outputs]


def test_customizer():
    atoms = bulk("Cu")
    results = bulk_to_slabs_flow(
        atoms, job_params={"static_job": {"asap_cutoff": True}}
    )
    for result in results:
        if result["name"] == "EMT Static":
            assert result["parameters"]["asap_cutoff"] is True


def test_customizer_v2():
    atoms = bulk("Cu")
    results = bulk_to_slabs_flow(atoms, job_params={"relax_job": {"asap_cutoff": True}})
    for result in results:
        if result["name"] == "EMT Static":
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
        if result["name"] == "EMT Static":
            assert result["parameters"]["asap_cutoff"] is False
        else:
            assert result["parameters"]["asap_cutoff"] is True
