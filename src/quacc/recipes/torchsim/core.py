from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, TypedDict

import torch_sim as ts
from torch_sim.autobatching import BinningAutoBatcher, InFlightAutoBatcher
from torch_sim.models.interface import ModelInterface

from quacc import job

if TYPE_CHECKING:
    import pathlib
    from collections.abc import Callable

    import torch
    from ase.atoms import Atoms
    from torch_sim.integrators import Integrator
    from torch_sim.optimizers import Optimizer
    from torch_sim.trajectory import TrajectoryReporter
    from torch_sim.typing import StateLike

    class TorchSimSchema(TypedDict):
        atoms: list[Atoms]
        model: str  # string of the model name

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
        max_steps: int
        steps_between_swaps: int
        trajectory_reporter: TrajectoryReporterDetails | None
        autobatcher: AutobatcherDetails
        init_kwargs: dict[str, Any] | None
        optimizer_kwargs: dict[str, Any]

    class TorchSimIntegrateSchema(TorchSimSchema):
        integrator: Integrator

        n_steps: int
        temperature: float | list[float]
        timestep: float
        trajectory_reporter: TrajectoryReporterDetails | None
        autobatcher: AutobatcherDetails
        integrator_kwargs: dict[str, Any]

    class TorchSimStaticSchema(TorchSimSchema):
        all_properties: list[dict[str, torch.Tensor]]
        trajectory_reporter: TrajectoryReporterDetails | None
        autobatcher: AutobatcherDetails


def _get_autobatcher_dict(
    autobatcher: InFlightAutoBatcher | BinningAutoBatcher,
) -> AutobatcherDetails:
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


def _get_reporter_dict(
    trajectory_reporter: TrajectoryReporter | dict | None,
    properties: list[str] | None = None,
) -> TrajectoryReporterDetails | None:
    if trajectory_reporter is not None:
        trajectory_reporter = ts.runners._configure_reporter(
            trajectory_reporter, properties=properties
        )

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
    else:
        return None


ModelType = Literal[
    "FairChemV1Model",
    "FairChemModel",
    "GraphPESWrapper",
    "MaceModel",
    "MatterSimModel",
    "MetatomicModel",
    "NequIPFrameworkModel",
    "OrbModel",
    "SevenNetModel",
]


def pick_model(model_type: ModelType, model, **model_kwargs) -> ModelInterface:
    if model_type == "FairChemV1Model":
        from torch_sim.models.fairchem_legacy import FairChemV1Model

        return FairChemV1Model(model=model, **model_kwargs)
    elif model_type == "FairChemModel":
        from torch_sim.models.fairchem import FairChemModel

        return FairChemModel(model=model, **model_kwargs)
    elif model_type == "GraphPESWrapper":
        from torch_sim.models.graphpes import GraphPESWrapper

        return GraphPESWrapper(model=model, **model_kwargs)
    elif model_type == "MaceModel":
        from torch_sim.models.mace import MaceModel

        return MaceModel(model=model, **model_kwargs)
    elif model_type == "MatterSimModel":
        from torch_sim.models.mattersim import MatterSimModel

        return MatterSimModel(model=model, **model_kwargs)
    elif model_type == "MetatomicModel":
        from torch_sim.models.metatomic import MetatomicModel

        return MetatomicModel(model=model, **model_kwargs)
    elif model_type == "NequIPFrameworkModel":
        from torch_sim.models.nequip_framework import NequIPFrameworkModel

        return NequIPFrameworkModel(model=model, **model_kwargs)
    elif model_type == "OrbModel":
        from torch_sim.models.orb import OrbModel

        return OrbModel(model=model, **model_kwargs)
    elif model_type == "SevenNetModel":
        from torch_sim.models.sevennet import SevenNetModel

        return SevenNetModel(model=model, **model_kwargs)
    else:
        raise ValueError(f"Invalid model type: {model_type}")


@job
def relax_job(
    atoms: list[Atoms],
    model: tuple[ModelType, Any] | ModelInterface,
    optimizer: Optimizer,
    convergence_fn: (
        Callable[[StateLike, torch.Tensor | None], torch.Tensor] | None
    ) = None,
    trajectory_reporter: TrajectoryReporter | dict | None = None,
    autobatcher: InFlightAutoBatcher | bool = False,
    max_steps: int = 10_000,
    steps_between_swaps: int = 5,
    init_kwargs: dict[str, Any] | None = None,
    **optimizer_kwargs: Any,
) -> TorchSimOptSchema:
    reporter_dict = _get_reporter_dict(
        trajectory_reporter, properties=["potential_energy"]
    )

    model = model if isinstance(model, ModelInterface) else pick_model(*model)

    state = ts.initialize_state(atoms, model.device, model.dtype)
    max_iterations = max_steps // steps_between_swaps
    autobatcher = ts.runners._configure_in_flight_autobatcher(
        state, model, autobatcher=autobatcher, max_iterations=max_iterations
    )
    autobatcher_dict = _get_autobatcher_dict(autobatcher)

    state = ts.optimize(
        system=state,
        model=model,
        optimizer=optimizer,
        convergence_fn=convergence_fn,
        trajectory_reporter=trajectory_reporter,
        autobatcher=autobatcher,
        max_steps=max_steps,
        steps_between_swaps=steps_between_swaps,
        init_kwargs=init_kwargs,
        **optimizer_kwargs,
    )

    return {
        "atoms": state.to_atoms(),
        "model": model.__class__.__name__,
        "optimizer": optimizer,
        "trajectory_reporter": reporter_dict,
        "autobatcher": autobatcher_dict,
        "max_steps": max_steps,
        "steps_between_swaps": steps_between_swaps,
        "init_kwargs": init_kwargs,
        "optimizer_kwargs": optimizer_kwargs,
    }


@job
def md_job(
    atoms: list[Atoms],
    model: tuple[ModelType, Any] | ModelInterface,
    integrator: Integrator,
    n_steps: int,
    temperature: float | list,
    timestep: float,
    trajectory_reporter: TrajectoryReporter | dict | None = None,
    autobatcher: BinningAutoBatcher | bool = False,
    **integrator_kwargs: Any,
) -> TorchSimIntegrateSchema:
    reporter_dict = _get_reporter_dict(
        trajectory_reporter,
        properties=["potential_energy", "kinetic_energy", "temperature"],
    )

    model = model if isinstance(model, ModelInterface) else pick_model(*model)

    state = ts.initialize_state(atoms, model.device, model.dtype)
    if autobatcher:
        autobatcher = ts.runners._configure_batches_iterator(
            state, model, autobatcher=autobatcher
        )  # type: ignore
        autobatcher_dict = _get_autobatcher_dict(autobatcher)
    else:
        autobatcher = False
        autobatcher_dict = None

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
        "model": model.__class__.__name__,
        "integrator": integrator,
        "n_steps": n_steps,
        "temperature": temperature,
        "timestep": timestep,
        "trajectory_reporter": reporter_dict,
        "autobatcher": autobatcher_dict,
        "integrator_kwargs": integrator_kwargs,
    }


@job
def static_job(
    atoms: list[Atoms],
    model: tuple[ModelType, Any] | ModelInterface,
    trajectory_reporter: TrajectoryReporter | dict | None = None,
    autobatcher: BinningAutoBatcher | bool = False,
) -> TorchSimStaticSchema:
    reporter_dict = _get_reporter_dict(
        trajectory_reporter, properties=["potential_energy"]
    )

    model = model if isinstance(model, ModelInterface) else pick_model(*model)

    state = ts.initialize_state(atoms, model.device, model.dtype)
    if autobatcher:
        autobatcher = ts.runners._configure_batches_iterator(
            state, model, autobatcher=autobatcher
        )  # type: ignore
        autobatcher_dict = _get_autobatcher_dict(autobatcher)
    else:
        autobatcher = False
        autobatcher_dict = None

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
        "model": model.__class__.__name__,
        "trajectory_reporter": reporter_dict,
        "autobatcher": autobatcher_dict,
    }
