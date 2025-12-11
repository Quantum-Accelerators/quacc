"""Core recipes for TorchSim."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING, Any, Literal

from monty.dev import requires

from quacc import __version__, job
from quacc.recipes.torchsim._base import (
    pick_model,
    process_binning_autobatcher_dict,
    process_in_flight_autobatcher_dict,
    process_trajectory_reporter_dict,
)
from quacc.runners._base import BaseRunner
from quacc.schemas.torchsim import CONVERGENCE_FN_REGISTRY, TSModelType

has_torchsim = bool(find_spec("torch_sim"))
if has_torchsim:
    import torch_sim as ts

if TYPE_CHECKING:
    from pathlib import Path

    from ase.atoms import Atoms

    if has_torchsim:
        from torch_sim.integrators import Integrator
        from torch_sim.optimizers import Optimizer

    from quacc.recipes.torchsim._base import (
        AutobatcherDict,
        TorchSimIntegrateSchema,
        TorchSimOptSchema,
        TorchSimStaticSchema,
        TrajectoryReporterDict,
    )


@job
@requires(has_torchsim, "torch_sim is required for this function")
def relax_job(
    atoms: list[Atoms],
    model_type: TSModelType,
    model_path: str | Path,
    optimizer: Optimizer,
    *,
    convergence_fn: Literal["energy", "force"] = "force",
    trajectory_reporter_dict: TrajectoryReporterDict | None = None,
    autobatcher_dict: AutobatcherDict | bool = False,
    max_steps: int = 10_000,
    steps_between_swaps: int = 5,
    init_kwargs: dict[str, Any] | None = None,
    model_kwargs: dict[str, Any] | None = None,
    convergence_fn_kwargs: dict[str, Any] | None = None,
    **optimizer_kwargs: Any,
) -> TorchSimOptSchema:
    """
    Carry out a geometry optimization on a set of atoms.

    Parameters
    ----------
    atoms : list[Atoms]
        The list of atoms objects.
    model_type : TSModelType
        The type of model to use, limited to the types supported by TorchSim.
    model_path : str | Path
        The path to the model file or checkpoint.
    optimizer : Optimizer
        The TorchSim optimizer to use.
    convergence_fn : Literal["energy", "force"]
        The convergence function, either "energy" or "force". This will use either the
        ts.generate_energy_convergence_fn or ts.generate_force_convergence_fn function
        to interally generate the convergence function. Arguments can be supplied via
        the convergence_fn_kwargs argument. Used to select convergence function
        generators from quacc.schemas.torchsim.CONVERGENCE_FN_REGISTRY.
    trajectory_reporter_dict : TrajectoryReporterDict | None
        This dictionary defines the trajectory reporting behavior. This is a
        quacc-specific dictionary that allows for the configuration of the
        TrajectoryReporter. For a list of available keys, refer to the TorchSim
        TrajectoryReporter documentation or the TrajectoryReporterDict type-hint.
    autobatcher_dict : AutobatcherDict | bool
        This dictionary defines the autobatcher behavior. This is a quacc-specific
        dictionary that allows for the configuration of the autobatcher. For a list of
        available keys, refer to the TorchSim Autobatcher's documentation or the
        AutobatcherDict type-hint.
    max_steps : int
        The maximum number of steps to run for each optimization.
    steps_between_swaps : int
        Number of steps to take before checking convergence and swapping out states.
    init_kwargs : dict[str, Any] | None
        Keyword arguments passed to the optimizer initialization function.
    model_kwargs : dict[str, Any] | None
        Keyword arguments passed to the model.
    convergence_fn_kwargs : dict[str, Any] | None
        Keyword arguments passed to the convergence function generator.
    **optimizer_kwargs : Any
        Keyword arguments to pass to the optimizer step function.

    Returns
    -------
    TorchSimOptSchema
        A dictionary representing the final atoms configuration and metadata geometry
        optimization schema. See the type-hint for the data structure.
    """
    runner = BaseRunner()
    runner.setup()

    model = pick_model(model_type, model_path, **model_kwargs or {})

    trajectory_reporter, trajectory_reporter_details = process_trajectory_reporter_dict(
        trajectory_reporter_dict, runner, n_systems=len(atoms)
    )

    max_iterations = max_steps // steps_between_swaps
    autobatcher, autobatcher_details = process_in_flight_autobatcher_dict(
        atoms, model, autobatcher_dict=autobatcher_dict, max_iterations=max_iterations
    )

    convergence_fn_obj = CONVERGENCE_FN_REGISTRY[convergence_fn](
        **convergence_fn_kwargs or {}
    )

    state = ts.optimize(
        system=atoms,
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
    runner.cleanup()

    return {
        "atoms": state.to_atoms(),
        "dir_name": str(runner.job_results_dir),
        "model_type": model_type,
        "model_path": model_path,
        "optimizer": optimizer,
        "convergence_fn": convergence_fn,
        "trajectory_reporter": trajectory_reporter_details,
        "autobatcher": autobatcher_details,
        "max_steps": max_steps,
        "steps_between_swaps": steps_between_swaps,
        "init_kwargs": init_kwargs,
        "model_kwargs": model_kwargs,
        "convergence_fn_kwargs": convergence_fn_kwargs,
        "optimizer_kwargs": optimizer_kwargs,
        "quacc_version": __version__,
    }


@job
@requires(has_torchsim, "torch_sim is required for this function")
def md_job(
    atoms: list[Atoms],
    model_type: TSModelType,
    model_path: str | Path,
    integrator: Integrator,
    *,
    n_steps: int,
    temperature: float | list,
    timestep: float,
    trajectory_reporter_dict: TrajectoryReporterDict | None = None,
    autobatcher_dict: AutobatcherDict | bool = False,
    model_kwargs: dict[str, Any] | None = None,
    **integrator_kwargs: Any,
) -> TorchSimIntegrateSchema:
    """
    Carry out a molecular dynamics calculation on a set of atoms.

    Parameters
    ----------
    atoms : list[Atoms]
        The list of atoms objects.
    model_type : TSModelType
        The type of model to use, limited to the types supported by TorchSim.
    model_path : str | Path
        The path to the model file or checkpoint.
    integrator : Integrator
        The TorchSim integrator to use.
    n_steps : int
        The maximum number of steps to run for each integration.
    temperature : float | list
        The temperature to use.
    timestep : float
        The timestep to use.
    trajectory_reporter_dict : TrajectoryReporterDict | None
        This dictionary defines the trajectory reporting behavior. This is a
        quacc-specific dictionary that allows for the configuration of the
        TrajectoryReporter. For a list of available keys, refer to the TorchSim
        TrajectoryReporter documentation or the TrajectoryReporterDict type-hint.
    autobatcher_dict : AutobatcherDict | bool
        This dictionary defines the autobatcher behavior. This is a quacc-specific
        dictionary that allows for the configuration of the autobatcher. For a list of
        available keys, refer to the TorchSim Autobatcher's documentation or the
        AutobatcherDict type-hint.
    model_kwargs : dict[str, Any] | None
        Keyword arguments passed to the model.
    **integrator_kwargs : Any
        Keyword arguments to pass to the integrator step function.

    Returns
    -------
    TorchSimIntegrateSchema
        A dictionary representing the final atoms configuration and metadata for the
        molecular dynamics schema. See the type-hint for the data structure.
    """
    runner = BaseRunner()
    runner.setup()

    model = pick_model(model_type, model_path, **model_kwargs or {})

    # Configure trajectory reporter
    trajectory_reporter, trajectory_reporter_details = process_trajectory_reporter_dict(
        trajectory_reporter_dict, runner, n_systems=len(atoms)
    )

    # Configure autobatcher
    autobatcher, autobatcher_details = process_binning_autobatcher_dict(
        atoms, model, autobatcher_dict=autobatcher_dict
    )

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
    runner.cleanup()

    return {
        "atoms": state.to_atoms(),
        "dir_name": str(runner.job_results_dir),
        "model_type": model_type,
        "model_path": model_path,
        "integrator": integrator,
        "n_steps": n_steps,
        "temperature": temperature,
        "timestep": timestep,
        "trajectory_reporter": trajectory_reporter_details,
        "autobatcher": autobatcher_details,
        "model_kwargs": model_kwargs,
        "integrator_kwargs": integrator_kwargs,
        "quacc_version": __version__,
    }


@job
@requires(has_torchsim, "torch_sim is required for this function")
def static_job(
    atoms: list[Atoms],
    model_type: TSModelType,
    model_path: str | Path,
    *,
    trajectory_reporter_dict: TrajectoryReporterDict | None = None,
    autobatcher_dict: AutobatcherDict | bool = False,
    model_kwargs: dict[str, Any] | None = None,
) -> TorchSimStaticSchema:
    """
    Carry out a static calculation on a set of atoms.

    Parameters
    ----------
    atoms : list[Atoms]
        The list of atoms objects.
    model_type : TSModelType
        The type of model to use, limited to the types supported by TorchSim.
    model_path : str | Path
        The path to the model file or checkpoint.
    trajectory_reporter_dict : TrajectoryReporterDict | None
        This dictionary defines the trajectory reporting behavior. This is a
        quacc-specific dictionary that allows for the configuration of the
        TrajectoryReporter. For a list of available keys, refer to the TorchSim
        TrajectoryReporter documentation or the TrajectoryReporterDict type-hint.
    autobatcher_dict : AutobatcherDict | bool
        This dictionary defines the autobatcher behavior. This is a quacc-specific
        dictionary that allows for the configuration of the autobatcher. For a list of
        available keys, refer to the TorchSim Autobatcher's documentation or the
        AutobatcherDict type-hint.
    model_kwargs : dict[str, Any] | None
        Keyword arguments passed to the model.

    Returns
    -------
    TorchSimStaticSchema
        A dictionary representing the final atoms configuration and metadata for the
        static calculation schema. See the type-hint for the data structure.
    """
    runner = BaseRunner()
    runner.setup()

    model = pick_model(model_type, model_path, **model_kwargs or {})

    trajectory_reporter, trajectory_reporter_details = process_trajectory_reporter_dict(
        trajectory_reporter_dict, runner, n_systems=len(atoms)
    )

    autobatcher, autobatcher_details = process_binning_autobatcher_dict(
        atoms, model, autobatcher_dict=autobatcher_dict
    )

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
    runner.cleanup()

    return {
        "atoms": atoms,
        "dir_name": str(runner.job_results_dir),
        "all_properties": all_properties_numpy,
        "model_type": model_type,
        "model_path": model_path,
        "trajectory_reporter": trajectory_reporter_details,
        "autobatcher": autobatcher_details,
        "model_kwargs": model_kwargs,
        "quacc_version": __version__,
    }
