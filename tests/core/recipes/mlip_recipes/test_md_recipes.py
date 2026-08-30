from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

import numpy as np
from ase.build import bulk
from ase.md.langevin import Langevin
from ase.units import fs

from quacc.recipes.mlip.md import md_job

libraries = []
if has_matcalc := find_spec("matcalc") and find_spec("matgl"):
    libraries.append("matcalc")


if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if get_token():
        libraries.append("fairchem")


def _calc_kwargs(library):
    if library == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        return {"name_or_path": "uma-s-1p1", "task_name": "omat"}
    if library == "matcalc":
        return {"name": "TensorNet-PES-MatPES-PBE-2025.2"}
    return {}


@pytest.mark.parametrize("library", libraries)
def test_md_job(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)
    old_positions = atoms.positions.copy()

    output = md_job(
        atoms,
        library=library,
        steps=20,
        md_params={
            "maxwell_boltzmann_kwargs": {
                "temperature_K": 300,
                "rng": np.random.default_rng(seed=42),
            },
            "set_com_stationary": True,
        },
        **_calc_kwargs(library),
    )

    assert output["name"] == f"{library} MD"
    assert len(output["trajectory"]) == 21
    assert output["parameters_md"]["timestep"] == pytest.approx(0.5 * fs)
    assert output["trajectory_log"][10]["time"] == pytest.approx(10 * 0.5 * fs)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    # The velocities are seeded, so the initial temperature does not depend on
    # the MLIP library.
    assert output["trajectory_log"][0]["temperature"] == pytest.approx(
        205.312, abs=0.01
    )
    # The input Atoms object should not be mutated.
    assert atoms.positions == pytest.approx(old_positions)


@pytest.mark.parametrize("library", libraries)
def test_md_job_nvt(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)

    output = md_job(
        atoms,
        library=library,
        dynamics=Langevin,
        steps=20,
        timestep_fs=0.5,
        temperature_K=300,
        md_params={
            "dynamics_kwargs": {"friction": 0.01},
            "maxwell_boltzmann_kwargs": {"rng": np.random.default_rng(seed=42)},
        },
        **_calc_kwargs(library),
    )

    assert output["name"] == f"{library} MD"
    assert len(output["trajectory"]) == 21
    assert output["parameters_md"]["timestep"] == pytest.approx(0.5 * fs)
    assert output["parameters_md"]["temperature_K"] == 300
    assert output["trajectory_log"][-1]["temperature"] > 0


@pytest.mark.parametrize("library", libraries)
def test_md_job_ensemble_nvt_langevin(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (2, 2, 2)

    output = md_job(
        atoms,
        library=library,
        dynamics="nvt_langevin",
        steps=10,
        temperature_K=300,
        md_params={"maxwell_boltzmann_kwargs": {"rng": np.random.default_rng(seed=42)}},
        **_calc_kwargs(library),
    )

    assert output["name"] == f"{library} MD"
    assert len(output["trajectory"]) == 11
    assert output["parameters_md"]["timestep"] == pytest.approx(0.5 * fs)
    assert output["parameters_md"]["temperature_K"] == 300
    assert output["parameters_md"]["friction"] == pytest.approx(1.0e-3)
    assert output["trajectory_log"][-1]["temperature"] > 0


@pytest.mark.parametrize("library", libraries)
def test_md_job_ensemble_npt(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    # The primitive fcc cell is not upper-triangular, which the ASE NPT class
    # requires; md_job should transform it automatically.
    atoms = bulk("Cu") * (2, 2, 2)
    old_cell = atoms.cell.copy()

    output = md_job(
        atoms,
        library=library,
        dynamics="npt",
        steps=10,
        temperature_K=300,
        md_params={"maxwell_boltzmann_kwargs": {"rng": np.random.default_rng(seed=42)}},
        **_calc_kwargs(library),
    )

    assert output["name"] == f"{library} MD"
    assert len(output["trajectory"]) == 11
    assert output["trajectory_log"][-1]["temperature"] > 0
    final_cell = output["atoms"].cell
    assert final_cell[1, 0] == final_cell[2, 0] == final_cell[2, 1] == 0.0
    # The input Atoms object should not be mutated.
    assert atoms.cell == pytest.approx(old_cell[:])
