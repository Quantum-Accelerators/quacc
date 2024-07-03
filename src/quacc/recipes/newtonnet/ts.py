"""Transition state recipes for the NewtonNet code."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires

from quacc import change_settings, get_settings, job, strip_decorator
from quacc.recipes.newtonnet.core import _add_stdev_and_hess, freq_job, relax_job
from quacc.runners.ase import Runner, run_neb
from quacc.schemas.ase import summarize_neb_run, summarize_opt_run
from quacc.utils.dicts import recursive_dict_merge

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))
has_sella = bool(find_spec("sella"))
has_newtonnet = bool(find_spec("newtonnet"))

if has_sella:
    from sella import IRC, Sella
if has_newtonnet:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
if has_geodesic_interpolate:
    from quacc.runners.ase import _geodesic_interpolate_wrapper

if TYPE_CHECKING:
    from typing import Any, Literal, TypedDict

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    from quacc.types import (
        NewtonNetIRCSchema,
        NewtonNetQuasiIRCSchema,
        NewtonNetTSSchema,
        OptParams,
        OptSchema,
    )

    class NebSchema(TypedDict):
        relax_reactant: OptSchema
        relax_product: OptSchema
        geodesic_results: list[Atoms]
        neb_results: dict

    class NebTsSchema(TypedDict):
        relax_reactant: OptSchema
        relax_product: OptSchema
        geodesic_results: list[Atoms]
        neb_results: dict
        ts_results: OptSchema


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def ts_job(
    atoms: Atoms,
    use_custom_hessian: bool = False,
    run_freq: bool = True,
    freq_job_kwargs: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    **calc_kwargs,
) -> NewtonNetTSSchema:
    """
    Perform a transition state (TS) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    run_freq
        Whether to run the frequency job.
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    TSSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    freq_job_kwargs = freq_job_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "optimizer": Sella,
        "optimizer_kwargs": (
            {"diag_every_n": 0, "order": 1} if use_custom_hessian else {"order": 1}
        ),
    }

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    if use_custom_hessian:
        opt_flags["optimizer_kwargs"]["hessian_function"] = _get_hessian

    calc = NewtonNet(**calc_flags)

    # Run the TS optimization
    dyn = Runner(atoms, calc).run_opt(**opt_flags)
    opt_ts_summary = _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"}), **calc_flags
    )

    # Run a frequency calculation
    freq_summary = (
        strip_decorator(freq_job)(opt_ts_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    opt_ts_summary["freq_job"] = freq_summary

    return opt_ts_summary


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    freq_job_kwargs: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    **calc_kwargs,
) -> NewtonNetIRCSchema:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    run_freq
        Whether to run the frequency analysis.
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    IRCSchema
        A dictionary containing the IRC summary and thermodynamic summary.
        See the type-hint for the data structure.
    """
    freq_job_kwargs = freq_job_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "optimizer": IRC,
        "optimizer_kwargs": {"dx": 0.1, "eta": 1e-4, "gamma": 0.4, "keep_going": True},
        "run_kwargs": {"direction": direction},
    }

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    # Define calculator
    calc = NewtonNet(**calc_flags)

    # Run IRC
    with change_settings({"CHECK_CONVERGENCE": False}):
        dyn = Runner(atoms, calc).run_opt(**opt_flags)
        opt_irc_summary = _add_stdev_and_hess(
            summarize_opt_run(
                dyn, additional_fields={"name": f"NewtonNet IRC: {direction}"}
            )
        )

    # Run frequency job
    freq_summary = (
        strip_decorator(freq_job)(opt_irc_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    opt_irc_summary["freq_job"] = freq_summary

    return opt_irc_summary


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def quasi_irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    irc_job_kwargs: dict[str, Any] | None = None,
    relax_job_kwargs: dict[str, Any] | None = None,
    freq_job_kwargs: dict[str, Any] | None = None,
) -> NewtonNetQuasiIRCSchema:
    """
    Perform a quasi-IRC job using the given atoms object. The initial IRC job by default
    is run with `max_steps: 5`.

    Parameters
    ----------
    atoms
        The atoms object representing the system
    direction
        The direction of the IRC calculation
    run_freq
        Whether to run the frequency analysis
    irc_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.irc_job][]
    relax_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.core.relax_job][]
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]

    Returns
    -------
    QuasiIRCSchema
        A dictionary containing the IRC summary, optimization summary, and
        thermodynamic summary.
        See the type-hint for the data structure.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    freq_job_kwargs = freq_job_kwargs or {}

    irc_job_defaults = {"max_steps": 5}
    irc_job_kwargs = recursive_dict_merge(irc_job_defaults, irc_job_kwargs)

    # Run IRC
    irc_summary = strip_decorator(irc_job)(
        atoms, direction=direction, run_freq=False, **irc_job_kwargs
    )

    # Run opt
    relax_summary = strip_decorator(relax_job)(irc_summary["atoms"], **relax_job_kwargs)

    # Run frequency
    freq_summary = (
        strip_decorator(freq_job)(relax_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    relax_summary["freq_job"] = freq_summary
    relax_summary["irc_job"] = irc_summary

    return relax_summary


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(
    has_geodesic_interpolate,
    "geodesic-interpolate must be installed. Refer to the quacc documentation.",
)
def neb_job(
    reactant_atoms: Atoms,
    product_atoms: Atoms,
    relax_job_kwargs: dict[str, Any] | None = None,
    calc_kwargs: dict[str, Any] | None = None,
    geodesic_interpolate_kwargs: dict[str, Any] | None = None,
    neb_kwargs: dict[str, Any] | None = None,
) -> NebSchema:
    """
    Perform a nudged elastic band (NEB) calculation to find the minimum energy path (MEP) between the given reactant and product structures.

    Parameters
    ----------
    reactant_atoms
        The Atoms object representing the reactant structure.
    product_atoms
        The Atoms object representing the product structure.
    relax_job_kwargs
        Keyword arguments to use relax_job.
    calc_kwargs
        Custom kwargs for the NewtonNet calculator.
    geodesic_interpolate_kwargs
        Keyword arguments for the geodesic function.
    neb_kwargs
        Keyword arguments for the NEB calculation.

    Returns
    -------
    NebSchema
        A dictionary containing the following keys:
            - 'relax_reactant': Summary of the relaxed reactant structure.
            - 'relax_product': Summary of the relaxed product structure.
            - 'geodesic_results': The interpolated images between reactant and product.
            - 'neb_results': Summary of the NEB optimization.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    neb_kwargs = neb_kwargs or {}
    geodesic_interpolate_kwargs = geodesic_interpolate_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": None,
    }

    geodesic_defaults = {"n_images": 20}

    neb_defaults = {"method": "aseneb", "precon": None}
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    geodesic_interpolate_flags = recursive_dict_merge(
        geodesic_defaults, geodesic_interpolate_kwargs
    )
    neb_flags = recursive_dict_merge(neb_defaults, neb_kwargs)

    # Debug prints to trace the values

    # Define calculator
    reactant_atoms.calc = NewtonNet(**calc_flags)
    product_atoms.calc = NewtonNet(**calc_flags)

    # Run relax job
    relax_summary_r = strip_decorator(relax_job)(reactant_atoms, **relax_job_kwargs)
    relax_summary_p = strip_decorator(relax_job)(product_atoms, **relax_job_kwargs)

    images = _geodesic_interpolate_wrapper(
        relax_summary_r["atoms"], relax_summary_p["atoms"], **geodesic_interpolate_flags
    )

    for image in images:
        image.calc = NewtonNet(**calc_flags)

    dyn = run_neb(images, neb_kwargs=neb_flags)

    return {
        "relax_reactant": relax_summary_r,
        "relax_product": relax_summary_p,
        "geodesic_results": images,
        "neb_results": summarize_neb_run(
            dyn,
            additional_fields={
                "neb_flags": neb_flags,
                "calc_flags": calc_flags,
                "geodesic_interpolate_flags": geodesic_interpolate_flags,
            },
        ),
    }


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(
    has_geodesic_interpolate,
    "geodesic-interpolate must be installed. Refer to the quacc documentation.",
)
def neb_ts_job(
    reactant_atoms: Atoms,
    product_atoms: Atoms,
    relax_job_kwargs: dict[str, Any] | None = None,
    calc_kwargs: dict[str, Any] | None = None,
    geodesic_interpolate_kwargs: dict[str, Any] | None = None,
    neb_kwargs: dict[str, Any] | None = None,
    ts_job_kwargs: dict[str, Any] | None = None,
) -> NebTsSchema:
    """
    Perform a quasi-IRC job using the given reactant and product atoms objects.

    Parameters
    ----------
    reactant_atoms
        The Atoms object representing the reactant structure.
    product_atoms
        The Atoms object representing the product structure.
    relax_job_kwargs
        Keyword arguments to use for the relax_job function, by default None.
    calc_kwargs
        Keyword arguments for the NewtonNet calculator, by default None.
    geodesic_interpolate_kwargs
        Keyword arguments for the geodesic_interpolate function, by default None.
    neb_kwargs
        Keyword arguments for the NEB calculation, by default None.
    ts_job_kwargs
        Keyword arguments for the TS optimizer, by default None.

    Returns
    -------
    NebTsSchema
        A dictionary containing the following keys:
            - 'relax_reactant': Summary of the relaxed reactant structure.
            - 'relax_product': Summary of the relaxed product structure.
            - 'geodesic_results': The interpolated images between reactant and product.
            - 'neb_results': Summary of the NEB optimization.
            - 'ts_results': Summary of the transition state optimization.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    neb_kwargs = neb_kwargs or {}
    geodesic_interpolate_kwargs = geodesic_interpolate_kwargs or {}
    calc_kwargs = calc_kwargs or {}
    ts_job_kwargs = ts_job_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": None,
    }

    geodesic_defaults = {"n_images": 20}

    neb_defaults = {"method": "aseneb", "precon": None}
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc_flags["hess_method"] = None
    geodesic_interpolate_flags = recursive_dict_merge(
        geodesic_defaults, geodesic_interpolate_kwargs
    )
    neb_flags = recursive_dict_merge(neb_defaults, neb_kwargs)

    neb_results = strip_decorator(neb_job)(
        reactant_atoms,
        product_atoms,
        calc_kwargs=calc_flags,
        geodesic_interpolate_kwargs=geodesic_interpolate_flags,
        neb_kwargs=neb_flags,
        relax_job_kwargs=relax_job_kwargs,
    )

    traj = neb_results["neb_results"]["trajectory"]
    traj_results = neb_results["neb_results"]["trajectory_results"]
    n_images = len(neb_results["geodesic_results"])

    ts_index = np.argmax([i["energy"] for i in traj_results[-(n_images - 1) : -1]]) + 1
    ts_atoms = traj[-(n_images) + ts_index]

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)

    output = strip_decorator(ts_job)(ts_atoms, **ts_job_kwargs, **calc_flags)
    neb_results["ts_results"] = output

    return neb_results


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(
    has_geodesic_interpolate,
    "geodesic-interpolate must be installed. Refer to the quacc documentation.",
)
def geodesic_job(
    reactant_atoms: Atoms,
    product_atoms: Atoms,
    relax_job_kwargs: dict[str, Any] | None = None,
    calc_kwargs: dict[str, Any] | None = None,
    geodesic_interpolate_kwargs: dict[str, Any] | None = None,
) -> dict:
    """
    Perform a quasi-IRC job using the given reactant and product atoms objects.

    Parameters
    ----------
    reactant_atoms
        The Atoms object representing the reactant structure.
    product_atoms
        The Atoms object representing the product structure.
    relax_job_kwargs
        Keyword arguments to use for the relax_job function, by default None.
    calc_kwargs
        Keyword arguments for the NewtonNet calculator, by default None.
    geodesic_interpolate_kwargs
        Keyword arguments for the geodesic_interpolate function, by default None.

    Returns
    -------
    dict
        A dictionary containing the following keys:
            - 'relax_reactant': Summary of the relaxed reactant structure.
            - 'relax_product': Summary of the relaxed product structure.
            - 'geodesic_results': The interpolated images between reactant and product.
            - 'highest_e_atoms': ASE atoms object for the highest energy structure for the geodesic path
    """
    relax_job_kwargs = relax_job_kwargs or {}
    geodesic_interpolate_kwargs = geodesic_interpolate_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": None,
    }

    geodesic_defaults = {"n_images": 20}

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc_flags["hess_method"] = None
    geodesic_interpolate_flags = recursive_dict_merge(
        geodesic_defaults, geodesic_interpolate_kwargs
    )

    # Define calculator
    reactant_atoms.calc = NewtonNet(**calc_flags)
    product_atoms.calc = NewtonNet(**calc_flags)

    # Run IRC
    relax_summary_r = strip_decorator(relax_job)(reactant_atoms, **relax_job_kwargs)
    relax_summary_p = strip_decorator(relax_job)(product_atoms, **relax_job_kwargs)

    images = _geodesic_interpolate_wrapper(
        relax_summary_r["atoms"].copy(),
        relax_summary_p["atoms"].copy(),
        **geodesic_interpolate_flags,
    )

    potential_energies = []
    for image in images:
        image.calc = NewtonNet(**calc_flags)
        potential_energies.append(image.get_potential_energy())

    ts_index = np.argmax(potential_energies)
    highest_e_atoms = images[ts_index]

    return {
        "relax_reactant": relax_summary_r,
        "relax_product": relax_summary_p,
        "geodesic_results": images,
        "highest_e_atoms": highest_e_atoms,
    }


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(
    has_geodesic_interpolate,
    "geodesic-interpolate must be installed. Refer to the quacc documentation.",
)
def geodesic_ts_job(
    reactant_atoms: Atoms,
    product_atoms: Atoms,
    relax_job_kwargs: dict[str, Any] | None = None,
    calc_kwargs: dict[str, Any] | None = None,
    geodesic_interpolate_kwargs: dict[str, Any] | None = None,
    ts_job_kwargs: dict[str, Any] | None = None,
) -> NebTsSchema:
    """
    Perform a quasi-IRC job using the given reactant and product atoms objects.

    Parameters
    ----------
    reactant_atoms
        The Atoms object representing the reactant structure.
    product_atoms
        The Atoms object representing the product structure.
    relax_job_kwargs
        Keyword arguments to use for the relax_job function, by default None.
    calc_kwargs
        Keyword arguments for the NewtonNet calculator, by default None.
    geodesic_interpolate_kwargs
        Keyword arguments for the geodesic_interpolate function, by default None.
    ts_job_kwargs
        Keyword arguments for ts optimizer, by default None.

    Returns
    -------
    NebTsSchema
        A dictionary containing the following keys:
            - 'relax_reactant': Summary of the relaxed reactant structure.
            - 'relax_product': Summary of the relaxed product structure.
            - 'geodesic_results': The interpolated images between reactant and product.
            - 'neb_results': Summary of the NEB optimization.
            - 'ts_results': Summary of the transition state optimization.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    ts_job_kwargs = ts_job_kwargs or {}
    geodesic_interpolate_kwargs = geodesic_interpolate_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": None,
    }

    geodesic_defaults = {"n_images": 20}

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc_flags["hess_method"] = None
    geodesic_interpolate_flags = recursive_dict_merge(
        geodesic_defaults, geodesic_interpolate_kwargs
    )

    # Define calculator
    reactant_atoms.calc = NewtonNet(**calc_flags)
    product_atoms.calc = NewtonNet(**calc_flags)

    # Run IRC
    relax_summary_r = strip_decorator(relax_job)(reactant_atoms, **relax_job_kwargs)
    relax_summary_p = strip_decorator(relax_job)(product_atoms, **relax_job_kwargs)

    images = _geodesic_interpolate_wrapper(
        relax_summary_r["atoms"].copy(),
        relax_summary_p["atoms"].copy(),
        **geodesic_interpolate_flags,
    )

    potential_energies = []
    for image in images:
        image.calc = NewtonNet(**calc_flags)
        potential_energies.append(image.get_potential_energy())

    ts_index = np.argmax(potential_energies)
    ts_atoms = images[ts_index]

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    output = strip_decorator(ts_job)(ts_atoms, **ts_job_kwargs, **calc_flags)
    return {
        "relax_reactant": relax_summary_r,
        "relax_product": relax_summary_p,
        "geodesic_results": images,
        "ts_results": output,
    }


def _get_hessian(atoms: Atoms) -> NDArray:
    """
    Calculate and retrieve the Hessian matrix for the given molecular configuration.

    This function takes an ASE Atoms object representing a molecular
    configuration and uses the NewtonNet machine learning calculator to
    calculate the Hessian matrix. The calculated Hessian matrix is then
    returned.

    Parameters
    ----------
    atoms
        The ASE Atoms object representing the molecular configuration.

    Returns
    -------
    NDArray
        The calculated Hessian matrix, reshaped into a 2D array.
    """
    settings = get_settings()
    ml_calculator = NewtonNet(
        model_path=settings.NEWTONNET_MODEL_PATH,
        settings_path=settings.NEWTONNET_CONFIG_PATH,
        hess_method="autograd",
    )
    ml_calculator.calculate(atoms)

    return ml_calculator.results["hessian"].reshape((-1, 3 * len(atoms)))
