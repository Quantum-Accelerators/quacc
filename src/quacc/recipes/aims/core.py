"""Core recipes for FHI-aims."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from ase.optimize import BFGS

from quacc import job
from quacc.atoms.core import check_is_metal
from quacc.recipes.aims._base import run_and_summarize, run_and_summarize_opt

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, RunSchema, SourceDirectory

BASE_SET_METAL = {
    "occupation_type": "cold 0.1",
    "relativistic": "atomic_zora scalar",
    "sc_accuracy_rho": 1e-5,
    "charge_mix_param": 0.05,
    "mixer": "pulay",
    "n_max_pulay": 14,
    "xc": "pbe",
    "output_level": "normal",
}

BASE_SET_NON_METAL = {
    "occupation_type": "gaussian 0.01",
    "relativistic": "atomic_zora scalar",
    "sc_accuracy_rho": 1e-5,
    "xc": "pbe",
    "charge_mix_param": 0.20,
    "mixer": "pulay",
    "output_level": "normal",
}

KSPACING_METAL = 0.033
KSPACING_NON_METAL = 0.045


@job
def static_job(
    atoms: Atoms,
    species_defaults: Literal[
        "light", "intermediate", "tight", "really_tight"
    ] = "intermediate",
    kspacing: float | None = None,
    spin: Literal["none", "collinear", "non-collinear"] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic SCF calculation with FHI-aims.

    Parameters
    ----------
    atoms
        The Atoms object.
    species_defaults
        The level of accuracy for the basis set and integration grids.
        Options: "light", "intermediate", "tight", "really_tight".
        Default is "intermediate" which is suitable for production SCF calculations.
    kspacing
        The kpoint spacing in Å^-1. If not provided, defaults to
        KSPACING_METAL (0.033) for metals and KSPACING_NON_METAL (0.045) for non-metals.
        Ignored for aperiodic systems.
    spin
        Spin treatment. Options are "none", "collinear", or "non-collinear".
        Default is None, which will automatically set to "collinear" if magnetic
        moments are detected in the atoms object.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the FHI-aims calculator. For a list of available
        keys, refer to the [ase.calculators.aims.Aims][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = (
        BASE_SET_METAL.copy() if check_is_metal(atoms) else BASE_SET_NON_METAL.copy()
    )
    calc_defaults["species_dir"] = species_defaults

    if kspacing is not None:
        calc_defaults["kspacing"] = kspacing
    else:
        calc_defaults["kspacing"] = (
            KSPACING_METAL if check_is_metal(atoms) else KSPACING_NON_METAL
        )

    if spin is not None:
        calc_defaults["spin"] = spin

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "FHI-aims Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    species_defaults: Literal[
        "light", "intermediate", "tight", "really_tight"
    ] = "light",
    kspacing: float | None = None,
    spin: Literal["none", "collinear", "non-collinear"] | None = None,
    relax_cell: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a structure relaxation with FHI-aims internal optimizer.

    Parameters
    ----------
    atoms
        The Atoms object.
    species_defaults
        The level of accuracy for the basis set and integration grids.
        Options: "light", "intermediate", "tight", "really_tight".
        Default is "light" which is suitable for geometry optimizations.
    kspacing
        The kpoint spacing in Å^-1. If not provided, defaults to
        KSPACING_METAL (0.033) for metals and KSPACING_NON_METAL (0.045) for non-metals.
        Ignored for aperiodic systems.
    spin
        Spin treatment. Options are "none", "collinear", or "non-collinear".
        Default is None, which will automatically set to "collinear" if magnetic
        moments are detected in the atoms object.
    relax_cell
        Whether to relax the cell or not.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the FHI-aims calculator. For a list of available
        keys, refer to the [ase.calculators.aims.Aims][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = (
        BASE_SET_METAL.copy() if check_is_metal(atoms) else BASE_SET_NON_METAL.copy()
    )
    calc_defaults["species_dir"] = species_defaults

    if kspacing is not None:
        calc_defaults["kspacing"] = kspacing
    else:
        calc_defaults["kspacing"] = (
            KSPACING_METAL if check_is_metal(atoms) else KSPACING_NON_METAL
        )

    if spin is not None:
        calc_defaults["spin"] = spin

    calc_defaults["relax_geometry"] = "bfgs 1E-2"

    if relax_cell:
        calc_defaults["relax_unit_cell"] = "full"

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "FHI-aims Relax"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def ase_relax_job(
    atoms: Atoms,
    species_defaults: Literal[
        "light", "intermediate", "tight", "really_tight"
    ] = "light",
    kspacing: float | None = None,
    spin: Literal["none", "collinear", "non-collinear"] | None = None,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a structure relaxation with FHI-aims using ASE
    external optimizers.

    Parameters
    ----------
    atoms
        The Atoms object.
    species_defaults
        The level of accuracy for the basis set and integration grids.
        Options: "light", "intermediate", "tight", "really_tight".
        Default is "light" which is suitable for geometry optimizations.
    kspacing
        The kpoint spacing in Å^-1. If not provided, defaults to
        KSPACING_METAL (0.033) for metals and KSPACING_NON_METAL (0.045) for non-metals.
        Ignored for aperiodic systems.
    spin
        Spin treatment. Options are "none", "collinear", or "non-collinear".
        Default is None, which will automatically set to "collinear" if magnetic
        moments are detected in the atoms object.
    relax_cell
        True if a volume relaxation should be performed. False if only the positions
        should be updated.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Additional keyword arguments to pass to the FHI-aims calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = (
        BASE_SET_METAL.copy() if check_is_metal(atoms) else BASE_SET_NON_METAL.copy()
    )
    calc_defaults["species_dir"] = species_defaults

    if kspacing is not None:
        calc_defaults["kspacing"] = kspacing
    else:
        calc_defaults["kspacing"] = (
            KSPACING_METAL if check_is_metal(atoms) else KSPACING_NON_METAL
        )

    if spin is not None:
        calc_defaults["spin"] = spin

    calc_defaults["compute_forces"] = True

    opt_defaults = {"optimizer": BFGS, "relax_cell": relax_cell}

    return run_and_summarize_opt(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "FHI-aims ASE Relax"} | (additional_fields or {}),
        copy_files=copy_files,
    )
