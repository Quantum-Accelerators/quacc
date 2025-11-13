from __future__ import annotations

import traceback

import pytest
from ase.atoms import Atoms
from ase.build import bulk

from quacc.schemas.torchsim import ConvergenceFn, TSModelType

try:
    import torch
    import torch_sim as ts
    from torch_sim.autobatching import InFlightAutoBatcher
    from torch_sim.models.lennard_jones import LennardJonesModel
    from torch_sim.models.mace import MaceModel, MaceUrls
except ImportError:
    pytest.skip(
        f"Torch-Sim not installed: {traceback.format_exc()}", allow_module_level=True
    )

from quacc.recipes.torchsim.core import md_job, relax_job, static_job

try:
    from mace.calculators.foundations_models import mace_mp
    from torch_sim.models.mace import MaceModel
except (ImportError, ValueError):
    pytest.skip(
        f"MACE not installed: {traceback.format_exc()}", allow_module_level=True
    )


@pytest.fixture
def raw_mace_mp():
    return mace_mp(model=MaceUrls.mace_mp_small, return_raw_model=True)


@pytest.fixture
def mace_model_path(raw_mace_mp):
    return raw_mace_mp.model_paths


@pytest.fixture
def ar_atoms() -> Atoms:
    """Create a face-centered cubic (FCC) Argon structure."""
    return bulk("Ar", "fcc", a=5.26, cubic=True)


@pytest.fixture
def fe_atoms() -> Atoms:
    """Create crystalline iron using ASE."""
    return bulk("Fe", "fcc", a=5.26, cubic=True)


@pytest.fixture
def ts_mace_model(raw_mace_mp) -> MaceModel:
    return MaceModel(
        model=raw_mace_mp,
        device=torch.device("cpu"),
        dtype=torch.float32,
        compute_forces=True,
        compute_stress=True,
    )


@pytest.fixture
def lj_model() -> LennardJonesModel:
    """Create a Lennard-Jones model with reasonable parameters for Ar."""
    return LennardJonesModel(
        use_neighbor_list=True,
        sigma=3.405,
        epsilon=0.0104,
        device=torch.device("cpu"),
        dtype=torch.float64,
        compute_forces=True,
        compute_stress=True,
        cutoff=2.5 * 3.405,
    )


def test_relax_job_comprehensive(ar_atoms: Atoms, tmp_path) -> None:
    """Test relax_job with all kwargs including trajectory reporter and autobatcher."""
    # Perturb the structure to make optimization meaningful
    ar_atoms.positions += 0.1

    n_systems = 2
    trajectory_reporter = {
        "filenames": [tmp_path / f"relax_{i}.h5md" for i in range(n_systems)],
        "state_frequency": 5,
        "prop_calculators": {1: {"potential_energy": lambda state: state.energy}},
    }

    # Create autobatcher
    autobatcher = False

    result = relax_job(
        atoms=[ar_atoms] * n_systems,
        model_type=TSModelType.LENNARD_JONES,
        model_path=None,
        optimizer=ts.Optimizer.fire,
        convergence_fn=ConvergenceFn.FORCE,
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
    assert all(f.is_file() for f in result["trajectory_reporter"]["filenames"])

    # Check autobatcher details
    assert result["autobatcher"] is None
    # assert result["autobatcher"]["autobatcher"] == "InFlightAutoBatcher"
    # assert result["autobatcher"]["memory_scales_with"] == "n_atoms"
    # assert result["autobatcher"]["max_memory_scaler"] == 8 + 1

    # Check other parameters
    assert result["max_steps"] == 500
    assert result["steps_between_swaps"] == 10
    assert result["init_kwargs"]["cell_filter"] == ts.CellFilter.unit


def test_relax_job_mace(
    ar_atoms: Atoms, lj_model: LennardJonesModel, tmp_path, raw_mace_mp
) -> None:
    """Test relax_job with all kwargs including trajectory reporter and autobatcher."""
    # Perturb the structure to make optimization meaningful
    ar_atoms.positions += 0.1

    n_systems = 2
    trajectory_reporter = {
        "filenames": [tmp_path / f"relax_{i}.h5md" for i in range(n_systems)],
        "state_frequency": 5,
        "prop_calculators": {1: {"potential_energy": lambda state: state.energy}},
    }

    # Create autobatcher
    autobatcher = InFlightAutoBatcher(
        model=lj_model, memory_scales_with="n_atoms", max_memory_scaler=260
    )

    # Create convergence function
    convergence_fn = ts.generate_force_convergence_fn(force_tol=1e-1)

    relax_job(
        atoms=[ar_atoms] * n_systems,
        model=("MaceModel", raw_mace_mp),
        optimizer=ts.Optimizer.fire,
        convergence_fn=convergence_fn,
        trajectory_reporter=trajectory_reporter,
        autobatcher=autobatcher,
        max_steps=500,
        steps_between_swaps=10,
        init_kwargs={"cell_filter": ts.CellFilter.unit},
    )


def test_md_job_comprehensive(ar_atoms: Atoms, tmp_path) -> None:
    """Test md_job with all kwargs including trajectory reporter and autobatcher."""
    n_systems = 2
    trajectory_reporter = {
        "filenames": [tmp_path / f"md_{i}.h5md" for i in range(n_systems)],
        "state_frequency": 2,
        "prop_calculators": {
            1: {
                "potential_energy": lambda state: state.energy,
                "kinetic_energy": lambda state: ts.calc_kinetic_energy(
                    momenta=state.momenta, masses=state.masses
                ),
                "temperature": lambda state: ts.calc_temperature(
                    momenta=state.momenta, masses=state.masses
                ),
            }
        },
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
    assert all(f.is_file() for f in result["trajectory_reporter"]["filenames"])

    # Check autobatcher details
    assert result["autobatcher"] is None
    # assert result["autobatcher"]["autobatcher"] == "BinningAutoBatcher"
    # assert result["autobatcher"]["memory_scales_with"] == "n_atoms"


def test_static_job_comprehensive(ar_atoms: Atoms, tmp_path) -> None:
    """Test static_job with all kwargs including trajectory reporter and autobatcher."""
    n_systems = 2
    trajectory_reporter = {
        "filenames": [tmp_path / f"static_{i}.h5md" for i in range(n_systems)],
        "state_frequency": 1,
        "prop_calculators": {1: {"potential_energy": lambda state: state.energy}},
        "state_kwargs": {"save_forces": True},
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
    assert result["trajectory_reporter"]["state_kwargs"]["save_forces"] is True
    assert all(f.is_file() for f in result["trajectory_reporter"]["filenames"])

    # Check autobatcher details
    assert result["autobatcher"] is None
    # assert result["autobatcher"]["autobatcher"] == "BinningAutoBatcher"
    # assert result["autobatcher"]["memory_scales_with"] == "n_atoms"
