from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

import numpy as np
from ase.build import bulk
from ase.md.andersen import Andersen
from ase.md.bussi import Bussi
from ase.md.langevin import Langevin
from ase.md.nose_hoover_chain import MTKNPT, IsotropicMTKNPT, NoseHooverChainNVT
from ase.md.npt import NPT
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen, NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.units import bar, fs

from quacc import Remove
from quacc.recipes.mlip.md import _resolve_ensemble, _upper_triangular_cell, md_job

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
    assert output["parameters_md"]["timestep"] == pytest.approx(1.0 * fs)
    assert output["trajectory_log"][10]["time"] == pytest.approx(10 * fs)
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
            "maxwell_boltzmann_kwargs": {
                "temperature_K": 300,
                "rng": np.random.default_rng(seed=42),
            },
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
        md_params={
            "maxwell_boltzmann_kwargs": {
                "temperature_K": 300,
                "rng": np.random.default_rng(seed=42),
            }
        },
        **_calc_kwargs(library),
    )

    assert output["name"] == f"{library} MD"
    assert len(output["trajectory"]) == 11
    assert output["parameters_md"]["timestep"] == pytest.approx(1.0 * fs)
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
        dynamics="npt_nose_hoover",
        steps=10,
        temperature_K=300,
        md_params={
            "maxwell_boltzmann_kwargs": {
                "temperature_K": 300,
                "rng": np.random.default_rng(seed=42),
            }
        },
        **_calc_kwargs(library),
    )

    assert output["name"] == f"{library} MD"
    assert len(output["trajectory"]) == 11
    assert output["trajectory_log"][-1]["temperature"] > 0
    final_cell = output["atoms"].cell
    assert final_cell[1, 0] == final_cell[2, 0] == final_cell[2, 1] == 0.0
    # The input Atoms object should not be mutated.
    assert atoms.cell == pytest.approx(old_cell[:])


@pytest.mark.parametrize(
    ("ensemble", "expected_class"),
    [
        ("nve", VelocityVerlet),
        ("nvt", NoseHooverChainNVT),
        ("nvt_nose_hoover", NoseHooverChainNVT),
        ("nvt_berendsen", NVTBerendsen),
        ("nvt_langevin", Langevin),
        ("nvt_andersen", Andersen),
        ("nvt_bussi", Bussi),
        ("npt", NPT),
        ("npt_nose_hoover", NPT),
        ("npt_berendsen", NPTBerendsen),
        ("npt_inhomogeneous", Inhomogeneous_NPTBerendsen),
        ("npt_mtk", MTKNPT),
        ("npt_isotropic_mtk", IsotropicMTKNPT),
    ],
)
def test_resolve_ensemble_classes(ensemble, expected_class):
    dynamics, _ = _resolve_ensemble(ensemble, 1.0 * fs, 300, None)
    assert dynamics is expected_class


def test_resolve_ensemble_is_case_insensitive():
    dynamics, _ = _resolve_ensemble("NVT_Langevin", 1.0 * fs, 300, None)
    assert dynamics is Langevin


def test_resolve_ensemble_coupling_times_scale_with_timestep():
    _, kwargs = _resolve_ensemble("npt_berendsen", 2.0 * fs, 300, 1.0)
    assert kwargs["taut"] == pytest.approx(200 * fs)
    assert kwargs["taup"] == pytest.approx(2000 * fs)

    _, kwargs = _resolve_ensemble("nvt", 2.0 * fs, 300, None)
    assert kwargs["tdamp"] == pytest.approx(200 * fs)


def test_resolve_ensemble_pressure_routing():
    # ase.md.npt.NPT takes `externalstress`; the Berendsen and MTK families
    # take `pressure_au`.
    _, kwargs = _resolve_ensemble("npt", 1.0 * fs, 300, 2.0)
    assert kwargs["externalstress"] == pytest.approx(2.0 * bar)
    assert "pressure_au" not in kwargs

    _, kwargs = _resolve_ensemble("npt_berendsen", 1.0 * fs, 300, 2.0)
    assert kwargs["pressure_au"] == pytest.approx(2.0 * bar)
    assert "externalstress" not in kwargs

    # The pressure defaults to 0 bar for the NPT-family ensembles.
    _, kwargs = _resolve_ensemble("npt_mtk", 1.0 * fs, 300, None)
    assert kwargs["pressure_au"] == 0.0


def test_resolve_ensemble_temperature_handling():
    # An explicit 0 K is honored; only None means unset.
    _, kwargs = _resolve_ensemble("nvt_langevin", 1.0 * fs, 0, None)
    assert kwargs["temperature_K"] == 0

    _, kwargs = _resolve_ensemble("nvt_langevin", 1.0 * fs, None, None)
    assert kwargs["temperature_K"] is Remove

    # NVE ignores the temperature and pressure entirely.
    _, kwargs = _resolve_ensemble("nve", 1.0 * fs, 300, 1.0)
    assert kwargs == {}


def test_resolve_ensemble_unknown():
    with pytest.raises(ValueError, match="Unsupported ensemble"):
        _resolve_ensemble("nvt_gibberish", 1.0 * fs, 300, None)


def test_upper_triangular_cell():
    # The primitive fcc cell is not upper-triangular.
    atoms = bulk("Cu")
    old_cell = atoms.cell.copy()

    transformed = _upper_triangular_cell(atoms)
    new_cell = transformed.cell

    assert transformed is not atoms
    assert new_cell[1, 0] == new_cell[2, 0] == new_cell[2, 1] == 0.0
    assert new_cell.cellpar() == pytest.approx(old_cell.cellpar())
    # The input Atoms object should not be mutated.
    assert atoms.cell == pytest.approx(old_cell[:])

    # An already upper-triangular cell is returned as-is.
    cubic = bulk("Cu", cubic=True)
    assert _upper_triangular_cell(cubic) is cubic
