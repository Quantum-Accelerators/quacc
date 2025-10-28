from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, TypedDict

import torch_sim as ts
from torch_sim.autobatching import BinningAutoBatcher, InFlightAutoBatcher

from quacc import job

if TYPE_CHECKING:
    import pathlib
    from collections.abc import Callable

    import torch
    from ase.atoms import Atoms
    from torch_sim.integrators import Integrator
    from torch_sim.models.interface import ModelInterface
    from torch_sim.optimizers import Optimizer
    from torch_sim.trajectory import TrajectoryReporter
    from torch_sim.typing import StateLike

    class TorchSimSchema(TypedDict):
        atoms: list[Atoms]
        model: str  # string of the model name

    class TrajectoryReporterDetails(TypedDict):
        state_frequency: int
        trajectory_kwargs: dict[str, Any]
        # prop_calculators: dict[int, dict[str, Callable]]
        prop_calculators: dict[int, list[str]]
        state_kwargs: dict[str, Any]
        metadata: dict[str, str] | None
        shape_warned: bool
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
        trajectory_reporter: TrajectoryReporterDetails | None
        autobatcher: AutobatcherDetails


def _get_autobatcher_dict(
    autobatcher: InFlightAutoBatcher | BinningAutoBatcher | bool,
    fallback_autobatcher: Literal["BinningAutoBatcher", "InFlightAutoBatcher"],
    fallback_memory_scales_with: Literal["n_atoms", "n_atoms_x_density"],
) -> AutobatcherDetails:
    if isinstance(autobatcher, InFlightAutoBatcher | BinningAutoBatcher):
        return {
            "autobatcher": type(autobatcher).__name__,  # type: ignore
            "memory_scales_with": autobatcher.memory_scales_with,  # type: ignore
            "max_memory_scaler": autobatcher.max_memory_scaler,
            "max_atoms_to_try": autobatcher.max_atoms_to_try,
            "memory_scaling_factor": autobatcher.memory_scaling_factor,
            "max_iterations": autobatcher.max_iterations,
            "max_memory_padding": autobatcher.max_memory_padding,
        }
    else:
        return {
            "autobatcher": fallback_autobatcher,
            "memory_scales_with": fallback_memory_scales_with,
            "max_memory_scaler": None,
            "max_atoms_to_try": None,
            "memory_scaling_factor": None,
            "max_iterations": None,
            "max_memory_padding": None,
        }


def _get_reporter_dict(
    trajectory_reporter: TrajectoryReporter | None,
) -> TrajectoryReporterDetails | None:
    if trajectory_reporter is not None:
        return {
            "state_frequency": trajectory_reporter.state_frequency,
            "trajectory_kwargs": trajectory_reporter.trajectory_kwargs,
            "prop_calculators": {
                i: list(calcs.keys())
                for i, calcs in trajectory_reporter.prop_calculators.items()
            },
            "state_kwargs": trajectory_reporter.state_kwargs,
            "metadata": trajectory_reporter.metadata,
            "shape_warned": trajectory_reporter.shape_warned,
            "filenames": trajectory_reporter.filenames,
        }
    else:
        return None


@job
def relax_job(
    atoms: Atoms,
    model: ModelInterface,
    optimizer: Optimizer,
    convergence_fn: (
        Callable[[StateLike, torch.Tensor | None], torch.Tensor] | None
    ) = None,
    trajectory_reporter: TrajectoryReporter | None = None,
    autobatcher: InFlightAutoBatcher | bool = False,
    max_steps: int = 10_000,
    steps_between_swaps: int = 5,
    init_kwargs: dict[str, Any] | None = None,
    **optimizer_kwargs: Any,
) -> TorchSimOptSchema:
    state = ts.optimize(
        system=atoms,
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

    reporter_dict = _get_reporter_dict(trajectory_reporter)
    autobatcher_dict = _get_autobatcher_dict(
        autobatcher, "InFlightAutoBatcher", model.memory_scales_with
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
    atoms: Atoms,
    model: ModelInterface,
    integrator: Integrator,
    n_steps: int,
    temperature: float | list,
    timestep: float,
    trajectory_reporter: TrajectoryReporter | None = None,
    autobatcher: BinningAutoBatcher | bool = False,
    **integrator_kwargs: Any,
) -> TorchSimIntegrateSchema:
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

    reporter_dict = _get_reporter_dict(trajectory_reporter)
    autobatcher_dict = _get_autobatcher_dict(
        autobatcher, "BinningAutoBatcher", model.memory_scales_with
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
    atoms: Atoms,
    model: ModelInterface,
    trajectory_reporter: TrajectoryReporter | None = None,
    autobatcher: BinningAutoBatcher | bool = False,
) -> TorchSimStaticSchema:
    state = ts.static(
        system=atoms,
        model=model,
        trajectory_reporter=trajectory_reporter,
        autobatcher=autobatcher,
    )

    reporter_dict = _get_reporter_dict(trajectory_reporter)
    autobatcher_dict = _get_autobatcher_dict(
        autobatcher, "BinningAutoBatcher", model.memory_scales_with
    )

    return {
        "atoms": state.to_atoms(),
        "model": model.__class__.__name__,
        "trajectory_reporter": reporter_dict,
        "autobatcher": autobatcher_dict,
    }
