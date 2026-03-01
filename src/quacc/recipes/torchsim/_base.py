"""Utility functions for TorchSim recipes."""

from __future__ import annotations

from copy import deepcopy
from importlib.util import find_spec
from typing import TYPE_CHECKING, Any, Literal, TypedDict

from monty.dev import requires

from quacc.schemas.torchsim import PROPERTY_FN_REGISTRY, TSModelType

has_torchsim = bool(find_spec("torch_sim"))

if has_torchsim:
    import torch_sim as ts
    from torch_sim.autobatching import BinningAutoBatcher, InFlightAutoBatcher


if TYPE_CHECKING:
    import pathlib
    from pathlib import Path

    from ase.atoms import Atoms

    if has_torchsim:
        from torch_sim.integrators import Integrator
        from torch_sim.models.interface import ModelInterface
        from torch_sim.optimizers import Optimizer
        from torch_sim.trajectory import TrajectoryReporter

    from quacc.runners._base import BaseRunner
    from quacc.schemas.torchsim import PropertyFn


class CalculationOutput(TypedDict):
    """Schema for the output of a TorchSim calculation."""

    energy: list[float]
    forces: list[list[list[float]]] | None
    stress: list[list[list[float]]] | None


class TorchSimSchema(TypedDict):
    atoms: list[Atoms]
    dir_name: str
    output: CalculationOutput
    model_type: TSModelType
    model_path: str | Path
    model_kwargs: dict[str, Any] | None
    trajectory_reporter: TrajectoryReporterDetails | None
    autobatcher: AutobatcherDetails | None
    quacc_version: str


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
    filenames: list[str | pathlib.Path] | None  # should be relative paths


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
    convergence_fn: Literal["energy", "force"]
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
    all_properties: list[dict[str, list]]


@requires(has_torchsim, "torch_sim is required for this function")
def process_in_flight_autobatcher_dict(
    atoms: list[Atoms],
    model: ModelInterface,
    autobatcher_dict: AutobatcherDict | bool,
    max_iterations: int,
) -> tuple[InFlightAutoBatcher | bool, AutobatcherDetails | None]:
    """Process the input dict into a InFlightAutoBatcher and details dictionary."""
    if isinstance(autobatcher_dict, bool):
        state = ts.initialize_state(atoms, model.device, model.dtype)
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


@requires(has_torchsim, "torch_sim is required for this function")
def process_binning_autobatcher_dict(
    atoms: list[Atoms], model: ModelInterface, autobatcher_dict: AutobatcherDict | bool
) -> tuple[BinningAutoBatcher | bool, AutobatcherDetails | None]:
    """Process the input dict into a BinningAutoBatcher and details dictionary."""
    if isinstance(autobatcher_dict, bool):
        # otherwise, configure the autobatcher, with the private runners method
        state = ts.initialize_state(atoms, model.device, model.dtype)
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


@requires(has_torchsim, "torch_sim is required for this function")
def process_trajectory_reporter_dict(
    trajectory_reporter_dict: TrajectoryReporterDict | None,
    runner: BaseRunner,
    n_systems: int,
) -> tuple[TrajectoryReporter, TrajectoryReporterDetails]:
    """Process the input dict into a TrajectoryReporter and details dictionary."""
    trajectory_reporter_dict = trajectory_reporter_dict or {}
    trajectory_reporter_dict = deepcopy(trajectory_reporter_dict)
    if "filenames" not in trajectory_reporter_dict:
        trajectory_reporter_dict["filenames"] = [
            runner.tmpdir / f"trajectory_{i}.h5md" for i in range(n_systems)
        ]
    else:
        trajectory_reporter_dict["filenames"] = [
            runner.tmpdir / filename
            for filename in trajectory_reporter_dict["filenames"]
        ]
    prop_calculators = trajectory_reporter_dict.pop("prop_calculators", {})
    prop_calculators_functions = {
        i: {prop: PROPERTY_FN_REGISTRY[prop] for prop in props}
        for i, props in prop_calculators.items()
    }

    # TODO: put in optional dependencies
    trajectory_reporter = ts.TrajectoryReporter(
        **trajectory_reporter_dict, prop_calculators=prop_calculators_functions
    )

    reporter_details = {
        "state_frequency": trajectory_reporter.state_frequency,
        "trajectory_kwargs": trajectory_reporter.trajectory_kwargs,
        "prop_calculators": prop_calculators,
        "state_kwargs": trajectory_reporter.state_kwargs,
        "metadata": trajectory_reporter.metadata,
        "filenames": trajectory_reporter.filenames,
    }
    return trajectory_reporter, reporter_details


@requires(has_torchsim, "torch_sim is required for this function")
def pick_model(
    model_type: TSModelType, model: str | Path, **model_kwargs: Any
) -> ModelInterface:
    """Pick and instantiate a model based on the model type.

    Parameters
    ----------
    model_type : TSModelType
        The type of model to instantiate.
    model : str | Path
        Path to the model file or checkpoint. For some models, string names may
        be allowed, such as "uma-s-1" for FairChemModel.
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

        return FairChemV1Model(model=model, **model_kwargs)
    if model_type == TSModelType.FAIRCHEM:
        from torch_sim.models.fairchem import FairChemModel

        return FairChemModel(model=model, **model_kwargs)
    if model_type == TSModelType.GRAPHPESWRAPPER:
        from torch_sim.models.graphpes import GraphPESWrapper

        return GraphPESWrapper(model=model, **model_kwargs)
    if model_type == TSModelType.MACE:
        from torch_sim.models.mace import MaceModel

        return MaceModel(model=model, **model_kwargs)
    if model_type == TSModelType.MATTERSIM:
        from torch_sim.models.mattersim import MatterSimModel

        return MatterSimModel(model=model, **model_kwargs)
    if model_type == TSModelType.METATOMIC:
        from torch_sim.models.metatomic import MetatomicModel

        return MetatomicModel(model=model, **model_kwargs)
    if model_type == TSModelType.NEQUIPFRAMEWORK:
        from torch_sim.models.nequip_framework import NequIPFrameworkModel

        return NequIPFrameworkModel(model=model, **model_kwargs)
    if model_type == TSModelType.ORB:
        from torch_sim.models.orb import OrbModel

        return OrbModel(model=model, **model_kwargs)
    if model_type == TSModelType.SEVENNET:
        from torch_sim.models.sevennet import SevenNetModel

        return SevenNetModel(model=model, **model_kwargs)
    if model_type == TSModelType.LENNARD_JONES:
        from torch_sim.models.lennard_jones import LennardJonesModel

        return LennardJonesModel(**model_kwargs)
    raise ValueError(f"Invalid model type: {model_type}")


def properties_to_calculation_output(
    all_properties_lists: list[dict[str, list]],
) -> CalculationOutput:
    """Convert properties from ts.static to a CalculationOutput.

    Parameters
    ----------
    all_properties_lists : list[dict[str, list]]
        List of property dictionaries from ts.static, with tensors converted to lists.

    Returns
    -------
    CalculationOutput
        The calculation output containing energy, forces, and stress.
    """
    energy = [prop_dict["potential_energy"][0] for prop_dict in all_properties_lists]
    forces = (
        [prop_dict["forces"] for prop_dict in all_properties_lists]
        if "forces" in all_properties_lists[-1]
        else None
    )
    stress = (
        [prop_dict["stress"][0] for prop_dict in all_properties_lists]
        if "stress" in all_properties_lists[-1]
        else None
    )
    return {"energy": energy, "forces": forces, "stress": stress}


@requires(has_torchsim, "torch_sim is required for this function")
def get_calculation_output(
    state: ts.SimState,
    model: ModelInterface,
    autobatcher: BinningAutoBatcher | InFlightAutoBatcher | bool = False,
) -> CalculationOutput:
    """Run a static calculation and return the output.

    Parameters
    ----------
    state : ts.SimState
        The simulation state to calculate properties for.
    model : ModelInterface
        The model to use for the calculation.
    autobatcher : BinningAutoBatcher | InFlightAutoBatcher | bool
        Optional autobatcher for batching calculations. If an InFlightAutoBatcher
        is passed, it will be converted to a BinningAutoBatcher.

    Returns
    -------
    CalculationOutput
        The calculation output containing energy, forces, and stress.
    """
    # Convert InFlightAutoBatcher to BinningAutoBatcher for ts.static
    if isinstance(autobatcher, InFlightAutoBatcher):
        autobatcher = BinningAutoBatcher(
            model=model,
            memory_scales_with=autobatcher.memory_scales_with,
            max_memory_scaler=autobatcher.max_memory_scaler,
        )

    properties = ts.static(system=state, model=model, autobatcher=autobatcher)

    all_properties_lists = [
        {name: t.tolist() for name, t in prop_dict.items()} for prop_dict in properties
    ]
    return properties_to_calculation_output(all_properties_lists)
