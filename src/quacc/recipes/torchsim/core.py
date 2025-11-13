from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, TypedDict

import torch_sim as ts
from torch_sim.autobatching import BinningAutoBatcher, InFlightAutoBatcher

from quacc import job
from quacc.schemas.torchsim import CONVERGENCE_FN_REGISTRY, ConvergenceFn, TSModelType

if TYPE_CHECKING:
    import pathlib
    from pathlib import Path

    import numpy as np
    from ase.atoms import Atoms
    from torch_sim.integrators import Integrator
    from torch_sim.models.interface import ModelInterface
    from torch_sim.optimizers import Optimizer
    from torch_sim.trajectory import TrajectoryReporter

    class TorchSimSchema(TypedDict):
        atoms: list[Atoms]
        model_type: TSModelType
        model_path: str | Path
        model_kwargs: dict[str, Any] | None
        trajectory_reporter: TrajectoryReporterDetails | None
        autobatcher: AutobatcherDetails | None

    class TrajectoryReporterDetails(TypedDict):
        state_frequency: int
        trajectory_kwargs: dict[str, Any]
        prop_calculators: dict[int, list[str]]
        state_kwargs: dict[str, Any]
        metadata: dict[str, str] | None
        filenames: list[str | pathlib.Path] | None

    class AutobatcherDetails(TypedDict):
        autobatcher: Literal["BinningAutoBatcher", "InFlightAutoBatcher"]
        memory_scales_with: Literal["n_atoms", "n_atoms_x_density"]
        max_memory_scaler: float | None
        max_atoms_to_try: int | None
        memory_scaling_factor: float | None
        max_iterations: int | None
        max_memory_padding: float | None

    class TorchSimOptSchema(TorchSimSchema):
        optimizer: Optimizer
        convergence_fn: ConvergenceFn
        max_steps: int
        steps_between_swaps: int
        init_kwargs: dict[str, Any] | None
        optimizer_kwargs: dict[str, Any]
        convergence_fn_kwargs: dict[str, Any] | None

    class TorchSimIntegrateSchema(TorchSimSchema):
        integrator: Integrator
        n_steps: int
        temperature: float | list[float]
        timestep: float
        integrator_kwargs: dict[str, Any]

    class TorchSimStaticSchema(TorchSimSchema):
        all_properties: list[dict[str, np.ndarray]]


def _get_autobatcher_details(
    autobatcher: InFlightAutoBatcher | BinningAutoBatcher,
) -> AutobatcherDetails:
    """Convert an autobatcher to a dictionary for serialization.

    Parameters
    ----------
    autobatcher : InFlightAutoBatcher | BinningAutoBatcher
        The autobatcher to convert.

    Returns
    -------
    AutobatcherDetails
        Dictionary representation of the autobatcher.
    """
    return {
        "autobatcher": type(autobatcher).__name__,  # type: ignore
        "memory_scales_with": autobatcher.memory_scales_with,  # type: ignore
        "max_memory_scaler": autobatcher.max_memory_scaler,
        "max_atoms_to_try": autobatcher.max_atoms_to_try,
        "memory_scaling_factor": autobatcher.memory_scaling_factor,
        "max_iterations": (
            autobatcher.max_iterations
            if isinstance(autobatcher, InFlightAutoBatcher)
            else None
        ),
        "max_memory_padding": autobatcher.max_memory_padding,
    }


def _get_reporter_details(
    trajectory_reporter: TrajectoryReporter | None,
) -> TrajectoryReporterDetails | None:
    """Convert a TrajectoryReporter to a dictionary for serialization.

    Parameters
    ----------
    trajectory_reporter : TrajectoryReporter | None
        The trajectory reporter to convert.

    Returns
    -------
    TrajectoryReporterDetails | None
        Dictionary representation of the trajectory reporter.
    """
    if trajectory_reporter is None:
        return None

    if trajectory_reporter.filenames is not None:
        filenames = [p.resolve() for p in trajectory_reporter.filenames]
    else:
        filenames = None

    return {
        "state_frequency": trajectory_reporter.state_frequency,
        "trajectory_kwargs": trajectory_reporter.trajectory_kwargs,
        "prop_calculators": {
            i: list(calcs.keys())
            for i, calcs in trajectory_reporter.prop_calculators.items()
        },
        "state_kwargs": trajectory_reporter.state_kwargs,
        "metadata": trajectory_reporter.metadata,
        "filenames": filenames,
    }


def pick_model(
    model_type: TSModelType, model_path: str | Path, **model_kwargs: Any
) -> ModelInterface:
    """Pick and instantiate a model based on the model type.

    Parameters
    ----------
    model_type : TSModelType
        The type of model to instantiate.
    model : str | Path
        Path to the model file or checkpoint.
    **model_kwargs : Any
        Additional keyword arguments to pass to the model constructor.

    Returns
    -------
    ModelInterface
        The instantiated model.

    Raises
    ------
    ValueError
        If an invalid model type is provided.
    """
    if model_type == TSModelType.FAIRCHEMV1:
        from torch_sim.models.fairchem_legacy import FairChemV1Model

        return FairChemV1Model(model=model_path, **model_kwargs)
    if model_type == TSModelType.FAIRCHEM:
        from torch_sim.models.fairchem import FairChemModel

        return FairChemModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.GRAPHPESWRAPPER:
        from torch_sim.models.graphpes import GraphPESWrapper

        return GraphPESWrapper(model=model_path, **model_kwargs)
    if model_type == TSModelType.MACE:
        from torch_sim.models.mace import MaceModel

        return MaceModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.MATTERSIM:
        from torch_sim.models.mattersim import MatterSimModel

        return MatterSimModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.METATOMIC:
        from torch_sim.models.metatomic import MetatomicModel

        return MetatomicModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.NEQUIPFRAMEWORK:
        from torch_sim.models.nequip_framework import NequIPFrameworkModel

        return NequIPFrameworkModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.ORB:
        from torch_sim.models.orb import OrbModel

        return OrbModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.SEVENNET:
        from torch_sim.models.sevennet import SevenNetModel

        return SevenNetModel(model=model_path, **model_kwargs)
    if model_type == TSModelType.LENNARD_JONES:
        from torch_sim.models.lennard_jones import LennardJonesModel

        return LennardJonesModel(**model_kwargs)
    raise ValueError(f"Invalid model type: {model_type}")


@job
def relax_job(
    atoms: list[Atoms],
    model_type: TSModelType,
    model_path: str | Path,
    optimizer: Optimizer,
    *,
    convergence_fn: ConvergenceFn = ConvergenceFn.FORCE,
    trajectory_reporter_dict: dict | None = None,
    autobatcher_dict: dict | bool = False,
    max_steps: int = 10_000,
    steps_between_swaps: int = 5,
    init_kwargs: dict[str, Any] | None = None,
    model_kwargs: dict[str, Any] | None = None,
    convergence_fn_kwargs: dict[str, Any] | None = None,
    **optimizer_kwargs: Any,
) -> TorchSimOptSchema:
    model = pick_model(model_type, model_path, **model_kwargs or {})

    convergence_fn_obj = CONVERGENCE_FN_REGISTRY[convergence_fn](
        **convergence_fn_kwargs or {}
    )

    state = ts.initialize_state(atoms, model.device, model.dtype)

    # Configure trajectory reporter
    trajectory_reporter = ts.runners._configure_reporter(
        trajectory_reporter_dict, properties=["potential_energy"]
    )

    # Configure autobatcher
    max_iterations = max_steps // steps_between_swaps
    if autobatcher_dict:
        autobatcher = ts.runners._configure_in_flight_autobatcher(
            state, model, autobatcher=autobatcher_dict, max_iterations=max_iterations
        )
    else:
        autobatcher = False

    state = ts.optimize(
        system=state,
        model=model,
        optimizer=optimizer,
        convergence_fn=convergence_fn_obj,
        trajectory_reporter=trajectory_reporter,
        autobatcher=autobatcher,
        max_steps=max_steps,
        steps_between_swaps=steps_between_swaps,
        init_kwargs=init_kwargs,
        **optimizer_kwargs,
    )

    return {
        "atoms": state.to_atoms(),
        "model_type": model_type,
        "model_path": model_path,
        "optimizer": optimizer,
        "convergence_fn": convergence_fn,
        "trajectory_reporter": _get_reporter_details(trajectory_reporter),
        "autobatcher": _get_autobatcher_details(autobatcher) if autobatcher else None,
        "max_steps": max_steps,
        "steps_between_swaps": steps_between_swaps,
        "init_kwargs": init_kwargs,
        "model_kwargs": model_kwargs,
        "convergence_fn_kwargs": convergence_fn_kwargs,
        "optimizer_kwargs": optimizer_kwargs,
    }


@job
def md_job(
    atoms: list[Atoms],
    model_type: TSModelType,
    model_path: str | Path,
    integrator: Integrator,
    *,
    n_steps: int,
    temperature: float | list,
    timestep: float,
    trajectory_reporter_dict: dict | None = None,
    autobatcher_dict: dict | bool = False,
    model_kwargs: dict[str, Any] | None = None,
    **integrator_kwargs: Any,
) -> TorchSimIntegrateSchema:
    model = pick_model(model_type, model_path, **model_kwargs or {})

    state = ts.initialize_state(atoms, model.device, model.dtype)

    # Configure trajectory reporter
    trajectory_reporter = ts.runners._configure_reporter(
        trajectory_reporter_dict,
        properties=["potential_energy", "kinetic_energy", "temperature"],
    )

    # Configure autobatcher
    if autobatcher_dict:
        autobatcher = ts.runners._configure_batches_iterator(
            state, model, autobatcher=autobatcher_dict
        )
    else:
        autobatcher = False

    state = ts.integrate(
        system=atoms,
        model=model,
        integrator=integrator,
        n_steps=n_steps,
        temperature=temperature,
        timestep=timestep,
        trajectory_reporter=trajectory_reporter,
        autobatcher=autobatcher,
        **integrator_kwargs,
    )

    return {
        "atoms": state.to_atoms(),
        "model_type": model_type,
        "model_path": model_path,
        "integrator": integrator,
        "n_steps": n_steps,
        "temperature": temperature,
        "timestep": timestep,
        "trajectory_reporter": _get_reporter_details(trajectory_reporter),
        "autobatcher": _get_autobatcher_details(autobatcher) if autobatcher else None,
        "model_kwargs": model_kwargs,
        "integrator_kwargs": integrator_kwargs,
    }


@job
def static_job(
    atoms: list[Atoms],
    model_type: TSModelType,
    model_path: str | Path,
    *,
    trajectory_reporter_dict: dict | None = None,
    autobatcher_dict: dict | bool = False,
    model_kwargs: dict[str, Any] | None = None,
) -> TorchSimStaticSchema:
    model = pick_model(model_type, model_path, **model_kwargs or {})

    state = ts.initialize_state(atoms, model.device, model.dtype)

    # Configure trajectory reporter
    trajectory_reporter = ts.runners._configure_reporter(
        trajectory_reporter_dict, properties=["potential_energy"]
    )

    # Configure autobatcher
    if autobatcher_dict:
        autobatcher = ts.runners._configure_batches_iterator(
            state, model, autobatcher=autobatcher_dict
        )
    else:
        autobatcher = False

    all_properties = ts.static(
        system=atoms,
        model=model,
        trajectory_reporter=trajectory_reporter,
        autobatcher=autobatcher,
    )

    all_properties_numpy = [
        {name: t.cpu().numpy() for name, t in prop_dict.items()}
        for prop_dict in all_properties
    ]

    return {
        "atoms": atoms,
        "all_properties": all_properties_numpy,
        "model_type": model_type,
        "model_path": model_path,
        "trajectory_reporter": _get_reporter_details(trajectory_reporter),
        "autobatcher": _get_autobatcher_details(autobatcher) if autobatcher else None,
        "model_kwargs": model_kwargs,
    }
