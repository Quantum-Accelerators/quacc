from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.optimize import FIRE

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.md import md_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True


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


def test_md_job1():
    atoms = molecule("H2O")
    old_positions = atoms.positions.copy()
    output = md_job(atoms, md_params={"steps": 500})
    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 501
    assert output["name"] == "EMT MD"
    assert output["parameters_md"]["timestep"] == pytest.approx(1.0)
    assert output["trajectory_log"][-1]["temperature"] == pytest.approx(1575.886)
    assert output["trajectory_log"][0]["temperature"] == pytest.approx(0.0)
    assert output["trajectory_log"][1]["temperature"] == pytest.approx(759.680)
    assert output["trajectory_log"][10]["time"] == pytest.approx(0.01)
    assert atoms.positions == pytest.approx(old_positions)


def test_md_job2():
    atoms = molecule("H2O")
    old_positions = atoms.positions.copy()

    output = md_job(
        atoms, md_params={"initial_temperature": 1000, "timestep": 0.5, "steps": 20}
    )
    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 21
    assert output["name"] == "EMT MD"
    assert output["parameters_md"]["timestep"] == pytest.approx(0.5)
    assert output["trajectory_log"][-1]["temperature"] == pytest.approx(1023.384)
    assert output["trajectory_log"][0]["temperature"] == pytest.approx(915.678)
    assert output["trajectory_log"][1]["temperature"] == pytest.approx(1060.650)
    assert output["trajectory_log"][10]["time"] == pytest.approx(0.005)
    assert atoms.positions == pytest.approx(old_positions)


def test_md_job3():
    atoms = molecule("H2O", vacuum=10.0)
    output = md_job(
        atoms,
        md_params={
            "timestep": 1.0,
            "dynamics": NPT,
            "dynamics_kwargs": {"temperature": 1000, "ttime": 50},
            "steps": 500,
        },
    )
    assert output["parameters"]["asap_cutoff"] is False
    assert len(output["trajectory"]) == 500
    assert output["name"] == "EMT MD"
    assert output["trajectory_log"][0]["temperature"] == pytest.approx(759.8829)
    assert output["trajectory_results"][-1]["energy"] == pytest.approx(2.0363759)


def test_md_job_error():
    atoms = molecule("H2O")
    with pytest.raises(ValueError, match="Quacc does not support"):
        md_job(atoms, md_params={"dynamics_kwargs": {"trajectory": "md.traj"}})


def test_md_logger(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O", vacuum=10.0)

    md_params = {
        "timestep": 1.0,
        "steps": 10,
        "dynamics": NPT,
        "dynamics_kwargs": {
            "temperature": 1000,
            "ttime": 50,
            "externalstress": 1,
            "pfactor": 40,
        },
    }

    with caplog.at_level(logging.WARNING):
        output_npt = md_job(atoms, md_params=md_params)

        assert (
            "In quacc `temperature`, `temp` and `temperature_K` are"
            " equivalent and will be interpreted in Kelvin."
        ) in caplog.text

    assert output_npt["parameters_md"]["md-type"] == "NPT"

    md_params = {
        "steps": 10,
        "dynamics": Langevin,
        "dynamics_kwargs": {
            "temperature_K": 1000,
            "friction": 0.01,
            "fixcm": False,
            "dt": 2.0,
        },
    }

    with caplog.at_level(logging.WARNING):
        output_langevin = md_job(atoms, md_params=md_params)

        assert "be interpreted as `timestep` in Quacc." in caplog.text

    assert output_langevin["parameters_md"]["md-type"] == "Langevin"
    assert output_langevin["parameters_md"]["timestep"] == pytest.approx(2.0)

    md_params = {
        "timestep": 1.0,
        "steps": 10,
        "dynamics": NVTBerendsen,
        "dynamics_kwargs": {"temperature": 1000, "taut": 10},
    }

    with caplog.at_level(logging.WARNING):
        output_nvt = md_job(atoms, md_params=md_params)

        assert (
            "In quacc `temperature`, `temp` and `temperature_K` are"
            " equivalent and will be interpreted in Kelvin."
        ) in caplog.text

    assert output_nvt["parameters_md"]["md-type"] == "NVTBerendsen"

    md_params = {
        "timestep": 1.0,
        "steps": 10,
        "dynamics": NPTBerendsen,
        "dynamics_kwargs": {
            "temperature": 1000,
            "pressure_au": 0.001,
            "taut": 500,
            "taup": 1000,
            "compressibility_au": 5e-4,
        },
    }

    with caplog.at_level(logging.WARNING):
        output_npt_berendsen = md_job(atoms, md_params=md_params)

        assert (
            "In quacc `temperature`, `temp` and `temperature_K` are"
            " equivalent and will be interpreted in Kelvin."
        ) in caplog.text

        assert (
            "In quacc `pressure` and `pressure_au` are equivalent"
            " and will be interpreted in GPA."
        ) in caplog.text

    assert output_npt_berendsen["parameters_md"]["md-type"] == "NPTBerendsen"

    md_params = {
        "timestep": 1.0,
        "steps": 10,
        "dynamics": NPTBerendsen,
        "dynamics_kwargs": {
            "temperature": 1000,
            "pressure": 0.001,
            "taut": 500,
            "taup": 1000,
            "compressibility": 5e-4,
        },
    }

    with caplog.at_level(logging.WARNING):
        output_npt_berendsen_bis = md_job(atoms, md_params=md_params)

        assert (
            "In quacc `temperature`, `temp` and `temperature_K` are"
            " equivalent and will be interpreted in Kelvin."
        ) in caplog.text

        assert (
            "In quacc `pressure` and `pressure_au` are equivalent"
            " and will be interpreted in GPA."
        ) in caplog.text

        assert (
            "In quacc `compressibility` and `compressibility_au` are equivalent"
            " and will be interpreted in 1/GPa."
        ) in caplog.text

    assert (
        output_npt_berendsen_bis["trajectory_log"]
        == output_npt_berendsen["trajectory_log"]
    )

    assert "ASE deprecated" in caplog.text
    assert "compressibility_au" in caplog.text
    assert "pressure_au" in caplog.text
    assert "temperature_K" in caplog.text
    assert "dt" in caplog.text


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
