from __future__ import annotations

import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import TYPE_CHECKING

import pytest
from ase.atoms import Atoms
from ase.build import bulk

torch = pytest.importorskip("torch")
ts = pytest.importorskip("torch_sim")
mace = pytest.importorskip("mace")

from mace.calculators.foundations_models import download_mace_mp_checkpoint

from quacc.recipes.torchsim import _base
from quacc.recipes.torchsim._base import (
    pick_model,
    process_binning_autobatcher_dict,
    process_in_flight_autobatcher_dict,
    process_trajectory_reporter_dict,
)
from quacc.recipes.torchsim.core import md_job, relax_job, static_job
from quacc.schemas.torchsim import TSModelType

if TYPE_CHECKING:
    if ts:
        from torchsim.core.models.lennard_jones_model import LennardJonesModel


@pytest.fixture
def mace_model_path():
    return Path(download_mace_mp_checkpoint("small"))


@pytest.fixture
def ar_atoms() -> Atoms:
    """Create a face-centered cubic (FCC) Argon structure."""
    return bulk("Ar", "fcc", a=5.26, cubic=True)


@pytest.fixture
def fe_atoms() -> Atoms:
    """Create crystalline iron using ASE."""
    return bulk("Fe", "fcc", a=5.26, cubic=True)


def test_autobatcher_helpers(monkeypatch, ar_atoms):
    class FakeAutoBatcher:
        memory_scales_with = "n_atoms"
        max_memory_scaler = 1
        max_atoms_to_try = 2
        memory_scaling_factor = 3
        max_iterations = 4
        max_memory_padding = 5

    model = SimpleNamespace(device="cpu", dtype="float64")
    monkeypatch.setattr(_base, "InFlightAutoBatcher", FakeAutoBatcher)
    monkeypatch.setattr(_base.ts, "initialize_state", lambda *args: object())
    monkeypatch.setattr(
        _base.ts.runners,
        "_configure_in_flight_autobatcher",
        lambda *args, **kwargs: FakeAutoBatcher(),
    )

    autobatcher, details = process_in_flight_autobatcher_dict(
        [ar_atoms], model, True, 4
    )
    assert isinstance(autobatcher, FakeAutoBatcher)
    assert details["max_iterations"] == 4

    monkeypatch.setattr(
        _base.ts.runners, "_configure_batches_iterator", lambda *args, **kwargs: []
    )
    assert process_binning_autobatcher_dict([ar_atoms], model, True) == (False, None)


def test_custom_trajectory_filenames(monkeypatch, tmp_path):
    class FakeReporter:
        def __init__(self, **kwargs):
            self.state_frequency = kwargs.get("state_frequency", 1)
            self.trajectory_kwargs = {}
            self.state_kwargs = {}
            self.metadata = None
            self.filenames = kwargs["filenames"]

    monkeypatch.setattr(_base.ts, "TrajectoryReporter", FakeReporter)
    runner = SimpleNamespace(tmpdir=tmp_path)

    _, details = process_trajectory_reporter_dict(
        {"filenames": ["custom.h5md"]}, runner, 1
    )

    assert details["filenames"] == [tmp_path / "custom.h5md"]


@pytest.mark.parametrize(
    ("model_type", "module_name", "class_name"),
    [
        (TSModelType.FAIRCHEMV1, "fairchem_legacy", "FairChemV1Model"),
        (TSModelType.FAIRCHEM, "fairchem", "FairChemModel"),
        (TSModelType.GRAPHPESWRAPPER, "graphpes", "GraphPESWrapper"),
        (TSModelType.MACE, "mace", "MaceModel"),
        (TSModelType.MATTERSIM, "mattersim", "MatterSimModel"),
        (TSModelType.METATOMIC, "metatomic", "MetatomicModel"),
        (TSModelType.NEQUIPFRAMEWORK, "nequip_framework", "NequIPFrameworkModel"),
        (TSModelType.ORB, "orb", "OrbModel"),
        (TSModelType.SEVENNET, "sevennet", "SevenNetModel"),
        (TSModelType.LENNARD_JONES, "lennard_jones", "LennardJonesModel"),
    ],
)
def test_pick_model(monkeypatch, model_type, module_name, class_name):
    class FakeModel:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    module = ModuleType(f"torch_sim.models.{module_name}")
    setattr(module, class_name, FakeModel)
    monkeypatch.setitem(sys.modules, module.__name__, module)

    model = pick_model(model_type, "model.pt", device="cpu")

    assert model.kwargs["device"] == "cpu"
    if model_type != TSModelType.LENNARD_JONES:
        assert model.kwargs["model"] == "model.pt"


def test_pick_model_invalid():
    with pytest.raises(ValueError, match="Invalid model type"):
        pick_model("invalid", "model.pt")


@pytest.fixture
def lj_model() -> LennardJonesModel:
    """Create a Lennard-Jones model with reasonable parameters for Ar."""
    return ts.LennardJonesModel(
        use_neighbor_list=True,
        sigma=3.405,
        epsilon=0.0104,
        device=torch.device("cpu"),
        dtype=torch.float64,
        compute_forces=True,
        compute_stress=True,
        cutoff=2.5 * 3.405,
    )


def test_relax_job_comprehensive(ar_atoms: Atoms) -> None:
    """Test relax_job with all kwargs including trajectory reporter and autobatcher."""
    # Perturb the structure to make optimization meaningful
    ar_atoms.positions += 0.1

    n_systems = 2
    trajectory_reporter = {
        "state_frequency": 5,
        "prop_calculators": {1: ["potential_energy"]},
    }

    # Create autobatcher
    autobatcher = False

    result = relax_job(
        atoms=[ar_atoms] * n_systems,
        model_type=TSModelType.LENNARD_JONES,
        model_path=None,
        optimizer=ts.Optimizer.fire,
        convergence_fn="force",
        trajectory_reporter_dict=trajectory_reporter,
        autobatcher_dict=autobatcher,
        max_steps=500,
        steps_between_swaps=10,
        init_kwargs={"cell_filter": ts.CellFilter.unit},
        model_kwargs={"sigma": 3.405, "epsilon": 0.0104, "compute_stress": True},
    )

    # Validate result structure
    assert "atoms" in result
    assert "model_type" in result
    assert "model_path" in result
    assert "optimizer" in result
    assert "convergence_fn" in result
    assert "trajectory_reporter" in result
    assert "autobatcher" in result
    assert "max_steps" in result
    assert "steps_between_swaps" in result
    assert "init_kwargs" in result
    assert "optimizer_kwargs" in result

    # Check atoms list output
    assert isinstance(result["atoms"], list)
    assert len(result["atoms"]) == n_systems
    assert isinstance(result["atoms"][0], Atoms)

    # Check model name
    assert result["model_type"] == TSModelType.LENNARD_JONES
    assert result["model_path"] is None

    # Check optimizer
    assert result["optimizer"] == ts.Optimizer.fire

    # Check trajectory reporter details
    assert result["trajectory_reporter"] is not None
    assert result["trajectory_reporter"]["state_frequency"] == 5
    assert "prop_calculators" in result["trajectory_reporter"]

    for i in range(n_systems):
        assert f"trajectory_{i}.h5md.gz" in os.listdir(result["dir_name"])

    # Check autobatcher details
    assert result["autobatcher"] is None

    # Check other parameters
    assert result["max_steps"] == 500
    assert result["steps_between_swaps"] == 10
    assert result["init_kwargs"]["cell_filter"] == ts.CellFilter.unit

    # Check calculation output
    assert "output" in result
    assert "energy" in result["output"]
    assert "forces" in result["output"]
    assert isinstance(result["output"]["energy"], list)
    assert len(result["output"]["energy"]) == n_systems
    assert isinstance(result["output"]["energy"][0], float)
    assert isinstance(result["output"]["forces"], list)
    assert len(result["output"]["forces"]) == n_systems


def test_relax_job_mace(ar_atoms: Atoms, mace_model_path: str) -> None:
    """Test relax_job with all kwargs including trajectory reporter and autobatcher."""
    # Perturb the structure to make optimization meaningful
    ar_atoms.positions += 0.1

    n_systems = 2
    trajectory_reporter_dict = {
        "state_frequency": 5,
        "prop_calculators": {1: ["potential_energy"]},
    }

    autobatcher_dict = {"memory_scales_with": "n_atoms", "max_memory_scaler": 260}

    relax_job(
        atoms=[ar_atoms] * n_systems,
        model_type=TSModelType.MACE,
        model_path=mace_model_path,
        optimizer=ts.Optimizer.fire,
        convergence_fn="force",
        trajectory_reporter_dict=trajectory_reporter_dict,
        autobatcher_dict=autobatcher_dict,
        max_steps=500,
        steps_between_swaps=10,
        init_kwargs={"cell_filter": ts.CellFilter.unit},
    )


def test_md_job_comprehensive(ar_atoms: Atoms) -> None:
    """Test md_job with all kwargs including trajectory reporter and autobatcher."""
    n_systems = 2
    trajectory_reporter = {
        "state_frequency": 2,
        "prop_calculators": {1: ["potential_energy", "kinetic_energy", "temperature"]},
    }

    # Create autobatcher
    autobatcher = False

    result = md_job(
        atoms=[ar_atoms] * n_systems,
        model_type=TSModelType.LENNARD_JONES,
        model_path=None,
        integrator=ts.Integrator.nvt_langevin,
        n_steps=20,
        temperature=300.0,
        timestep=0.001,
        trajectory_reporter_dict=trajectory_reporter,
        autobatcher_dict=autobatcher,
        model_kwargs={"sigma": 3.405, "epsilon": 0.0104, "compute_stress": True},
    )

    # Validate result structure
    assert "atoms" in result
    assert "model_type" in result
    assert "model_path" in result
    assert "integrator" in result
    assert "n_steps" in result
    assert "temperature" in result
    assert "timestep" in result
    assert "trajectory_reporter" in result
    assert "autobatcher" in result
    assert "integrator_kwargs" in result

    # Check atoms list output
    assert isinstance(result["atoms"], list)
    assert len(result["atoms"]) == n_systems
    assert isinstance(result["atoms"][0], Atoms)

    # Check model name
    assert result["model_type"] == TSModelType.LENNARD_JONES
    assert result["model_path"] is None

    # Check integrator
    assert result["integrator"] == ts.Integrator.nvt_langevin

    # Check MD parameters
    assert result["n_steps"] == 20
    assert result["temperature"] == 300.0
    assert result["timestep"] == 0.001

    # Check trajectory reporter details
    assert result["trajectory_reporter"] is not None
    assert result["trajectory_reporter"]["state_frequency"] == 2
    assert "prop_calculators" in result["trajectory_reporter"]
    for i in range(n_systems):
        assert f"trajectory_{i}.h5md.gz" in os.listdir(result["dir_name"])

    # Check autobatcher details
    assert result["autobatcher"] is None

    # Check calculation output
    assert "output" in result
    assert "energy" in result["output"]
    assert "forces" in result["output"]
    assert isinstance(result["output"]["energy"], list)
    assert len(result["output"]["energy"]) == n_systems
    assert isinstance(result["output"]["energy"][0], float)
    assert isinstance(result["output"]["forces"], list)
    assert len(result["output"]["forces"]) == n_systems


def test_static_job_comprehensive(ar_atoms: Atoms) -> None:
    """Test static_job with all kwargs including trajectory reporter and autobatcher."""
    n_systems = 2
    trajectory_reporter = {
        "state_frequency": 1,
        "prop_calculators": {1: ["potential_energy", "forces"]},
        "state_kwargs": {"save_forces": False},
    }

    # Create autobatcher
    autobatcher = False

    result = static_job(
        atoms=[ar_atoms] * n_systems,
        model_type=TSModelType.LENNARD_JONES,
        model_path=None,
        trajectory_reporter_dict=trajectory_reporter,
        autobatcher_dict=autobatcher,
        model_kwargs={"sigma": 3.405, "epsilon": 0.0104, "compute_stress": True},
    )

    # Validate result structure
    assert "atoms" in result
    assert "model_type" in result
    assert "model_path" in result
    assert "trajectory_reporter" in result
    assert "autobatcher" in result

    # Check atoms list output
    assert isinstance(result["atoms"], list)
    assert len(result["atoms"]) == n_systems
    assert isinstance(result["atoms"][0], Atoms)

    # Check model name
    assert result["model_type"] == TSModelType.LENNARD_JONES
    assert result["model_path"] is None

    # Check trajectory reporter details
    assert result["trajectory_reporter"] is not None
    assert result["trajectory_reporter"]["state_frequency"] == 1
    assert "prop_calculators" in result["trajectory_reporter"]
    assert "state_kwargs" in result["trajectory_reporter"]
    assert result["trajectory_reporter"]["state_kwargs"]["save_forces"] is False
    for i in range(n_systems):
        assert f"trajectory_{i}.h5md.gz" in os.listdir(result["dir_name"])

    # Check autobatcher details
    assert result["autobatcher"] is None

    # Check calculation output
    assert "output" in result
    assert "energy" in result["output"]
    assert "forces" in result["output"]
    assert isinstance(result["output"]["energy"], list)
    assert len(result["output"]["energy"]) == n_systems
    assert isinstance(result["output"]["energy"][0], float)
    assert isinstance(result["output"]["forces"], list)
    assert len(result["output"]["forces"]) == n_systems
