"""Utility functions for TorchSim recipes."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, TypedDict

import torch_sim as ts
from torch_sim.autobatching import BinningAutoBatcher, InFlightAutoBatcher

from quacc.schemas.torchsim import PROPERTY_FN_REGISTRY, ConvergenceFn, TSModelType

if TYPE_CHECKING:
    import pathlib

    import numpy as np
    from ase.atoms import Atoms
    from torch_sim.integrators import Integrator
    from torch_sim.models.interface import ModelInterface
    from torch_sim.optimizers import Optimizer
    from torch_sim.state import SimState
    from torch_sim.trajectory import TrajectoryReporter

    from quacc.schemas.torchsim import PropertyFn


class TorchSimSchema(TypedDict):
    atoms: list[Atoms]
    model_type: TSModelType
    model_path: str | Path
    model_kwargs: dict[str, Any] | None
    trajectory_reporter: TrajectoryReporterDetails | None
    autobatcher: AutobatcherDetails | None


class TrajectoryReporterDict(TypedDict, total=False):
    filenames: str | pathlib.Path | list[str | pathlib.Path]
    state_frequency: int | None
    prop_calculators: dict[int, list[PropertyFn]] | None
    state_kwargs: dict[str, Any] | None
    metadata: dict[str, str] | None
    trajectory_kwargs: dict[str, Any] | None


class AutobatcherDict(TypedDict, total=False):
    memory_scales_with: Literal["n_atoms", "n_atoms_x_density"] | None
    max_memory_scaler: float | None
    max_atoms_to_try: int | None
    memory_scaling_factor: float | None
    max_iterations: int | None
    max_memory_padding: float | None


class TrajectoryReporterDetails(TypedDict):
    state_frequency: int
    trajectory_kwargs: dict[str, Any]
    prop_calculators: dict[int, list[PropertyFn]]
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


def process_in_flight_autobatcher_dict(
    state: SimState,
    model: ModelInterface,
    autobatcher_dict: AutobatcherDict | bool,
    max_iterations: int,
) -> tuple[InFlightAutoBatcher | bool, AutobatcherDetails | None]:
    """Process the input dict into a InFlightAutoBatcher and details dictionary."""
    if isinstance(autobatcher_dict, bool):
        # False means no autobatcher
        if not autobatcher_dict:
            return False, None
        # otherwise, configure the autobatcher, with the private runners method
        autobatcher = ts.runners._configure_in_flight_autobatcher(
            state, model, autobatcher=autobatcher_dict, max_iterations=max_iterations
        )
    else:
        autobatcher = InFlightAutoBatcher(model=model, **autobatcher_dict)

    autobatcher_details = _get_autobatcher_details(autobatcher)
    return autobatcher, autobatcher_details


def process_binning_autobatcher_dict(
    state: SimState, model: ModelInterface, autobatcher_dict: AutobatcherDict | bool
) -> tuple[BinningAutoBatcher | bool, AutobatcherDetails | None]:
    """Process the input dict into a BinningAutoBatcher and details dictionary."""
    if isinstance(autobatcher_dict, bool):
        # otherwise, configure the autobatcher, with the private runners method
        autobatcher = ts.runners._configure_batches_iterator(
            state, model, autobatcher=autobatcher_dict
        )
        # list means no autobatcher
        if isinstance(autobatcher, list):
            return False, None
    else:
        # pop max_iterations if present
        autobatcher_dict.pop("max_iterations", None)
        autobatcher = BinningAutoBatcher(model=model, **autobatcher_dict)

    autobatcher_details = _get_autobatcher_details(autobatcher)
    return autobatcher, autobatcher_details


def _get_autobatcher_details(
    autobatcher: InFlightAutoBatcher | BinningAutoBatcher,
) -> AutobatcherDetails:
    """Extract the metadata of an autobatcher."""
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


def process_trajectory_reporter_dict(
    trajectory_reporter_dict: TrajectoryReporterDict | None,
) -> tuple[TrajectoryReporter, TrajectoryReporterDetails]:
    """Process the input dict into a TrajectoryReporter and details dictionary."""
    if trajectory_reporter_dict is None:
        return None, None
    trajectory_reporter_dict = deepcopy(trajectory_reporter_dict)

    prop_calculators = trajectory_reporter_dict.pop("prop_calculators", {})
    prop_calculators_functions = {
        i: {prop: PROPERTY_FN_REGISTRY[prop] for prop in props}
        for i, props in prop_calculators.items()
    }

    # TODO: put in optional dependencies
    trajectory_reporter = ts.TrajectoryReporter(
        **trajectory_reporter_dict, prop_calculators=prop_calculators_functions
    )

    trajectory_reporter.filenames = [
        Path(p).resolve() for p in trajectory_reporter_dict["filenames"]
    ]

    reporter_details = {
        "state_frequency": trajectory_reporter.state_frequency,
        "trajectory_kwargs": trajectory_reporter.trajectory_kwargs,
        "prop_calculators": prop_calculators,
        "state_kwargs": trajectory_reporter.state_kwargs,
        "metadata": trajectory_reporter.metadata,
        "filenames": trajectory_reporter.filenames,
    }
    return trajectory_reporter, reporter_details


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
