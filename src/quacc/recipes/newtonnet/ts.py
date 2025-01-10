"""Transition state recipes for the NewtonNet code."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
from ase.mep import NEB
from monty.dev import requires

from quacc import change_settings, get_settings, job, strip_decorator
from quacc.recipes.newtonnet.core import _add_stdev_and_hess, freq_job, relax_job
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))
has_sella = bool(find_spec("sella"))
has_newtonnet = bool(find_spec("newtonnet"))

if has_sella:
    from sella import IRC, Sella
if has_newtonnet:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
if has_geodesic_interpolate:
    from quacc.atoms.ts import geodesic_interpolate_wrapper

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    from quacc.types import (
        GeodesicSchema,
        NebSchema,
        NewtonNetIRCSchema,
        NewtonNetQuasiIRCSchema,
        NewtonNetTSSchema,
        OptParams,
    )


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
    additional_fields: dict[str, Any] | None = None,
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
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    TSSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    additional_fields = additional_fields or {}
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
        Summarize(additional_fields={"name": "NewtonNet TS"}).opt(dyn), **calc_flags
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
    additional_fields: dict[str, Any] | None = None,
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
    additional_fields
        Additional fields to add to the results dictionary.
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
            Summarize(
                additional_fields={"name": f"NewtonNet IRC: {direction}"}
                | (additional_fields or {})
            ).opt(dyn)
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
def neb_job(
    reactant_atoms: Atoms,
    product_atoms: Atoms,
    interpolation_method: Literal["linear", "idpp", "geodesic"] = "linear",
    relax_job_kwargs: dict[str, Any] | None = None,
    interpolate_kwargs: dict[str, Any] | None = None,
    neb_kwargs: dict[str, Any] | None = None,
    **calc_kwargs,
) -> NebSchema:
    """
    Perform a nudged elastic band (NEB) calculation to find the minimum energy path (MEP) between the given reactant and product structures.

    Parameters
    ----------
    reactant_atoms
        The Atoms object representing the reactant structure.
    product_atoms
        The Atoms object representing the product structure.
    interpolation_method
        The method to initialize the NEB optimization. There are three choices here, "linear", "idpp" and "geodesic".
        Defaults to linear.
    relax_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.core.relax_job][] function.
    interpolate_kwargs
        Keyword arguments for the interpolate functions (geodesic, linear or idpp).
    neb_kwargs
        Keyword arguments for the NEB calculation.
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    NebSchema
        A dictionary containing the NEB results
    """
    relax_job_kwargs = relax_job_kwargs or {}
    settings = get_settings()

    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": None,
    }

    interpolate_defaults = {"n_images": 20}

    neb_defaults = {"method": "aseneb", "precon": None}
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    interpolate_flags = recursive_dict_merge(interpolate_defaults, interpolate_kwargs)
    neb_flags = recursive_dict_merge(neb_defaults, neb_kwargs)

    # Define calculator
    calc = NewtonNet(**calc_flags)

    # Run relax job
    relax_summary_r = strip_decorator(relax_job)(reactant_atoms, **relax_job_kwargs)
    relax_summary_p = strip_decorator(relax_job)(product_atoms, **relax_job_kwargs)

    if interpolation_method == "geodesic":
        images = geodesic_interpolate_wrapper(
            relax_summary_r["atoms"], relax_summary_p["atoms"], **interpolate_flags
        )
    else:
        images = [reactant_atoms]
        images += [
            reactant_atoms.copy() for _ in range(interpolate_flags["n_images"] - 2)
        ]
        images += [product_atoms]
        neb = NEB(images)

        # Interpolate linearly the positions of the middle images:
        neb.interpolate(method=interpolation_method)
        images = neb.images

    max_steps = neb_flags.pop("max_steps", None)
    dyn = Runner(images, calc).run_neb(max_steps=max_steps, neb_kwargs=neb_flags)

    return {
        "relax_reactant": relax_summary_r,
        "relax_product": relax_summary_p,
        "initial_images": images,
        "neb_results": Summarize(
            additional_fields={
                "neb_flags": neb_flags,
                "calc_flags": calc_flags,
                "interpolate_flags": interpolate_flags,
            }
        ).neb(dyn, len(images)),
    }


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
    geodesic_interpolate_kwargs: dict[str, Any] | None = None,
    **calc_kwargs,
) -> GeodesicSchema:
    """
    Perform a quasi-IRC job using the given reactant and product atoms objects.

    Parameters
    ----------
    reactant_atoms
        The Atoms object representing the reactant structure.
    product_atoms
        The Atoms object representing the product structure.
    relax_job_kwargs
        Keyword arguments to use for [quacc.recipes.newtonnet.core.relax_job][] function.
    geodesic_interpolate_kwargs
        Keyword arguments for geodesic_interpolate, by default None.
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    GeodesicSchema
        A dictionary containing the following keys:
            - 'relax_reactant': Summary of the relaxed reactant structure.
            - 'relax_product': Summary of the relaxed product structure.
            - 'initial_images': The interpolated images between reactant and product.
            - 'ts_atoms': ASE atoms object for the highest energy structure for the geodesic path
    """
    relax_job_kwargs = relax_job_kwargs or {}
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

    images = geodesic_interpolate_wrapper(
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

    return {
        "relax_reactant": relax_summary_r,
        "relax_product": relax_summary_p,
        "initial_images": images,
        "ts_atoms": ts_atoms,
    }


def _get_hessian(atoms: Atoms, **calc_kwargs) -> NDArray:
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
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    NDArray
        The calculated Hessian matrix, reshaped into a 2D array.
    """
    settings = get_settings()
    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": "autograd",
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc = NewtonNet(**calc_flags)
    calc.calculate(atoms)
    return calc.results["hessian"].reshape((-1, 3 * len(atoms)))
