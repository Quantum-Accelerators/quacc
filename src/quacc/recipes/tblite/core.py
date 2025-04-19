"""Core recipes for the tblite code."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from ase.optimize import BFGS
from monty.dev import requires

from quacc import job
from quacc.recipes.common.core import Recipe
from quacc.schemas.thermo import ThermoSummarize
from quacc.utils.dicts import recursive_dict_merge

has_tblite = bool(find_spec("tblite"))
if has_tblite:
    from tblite.ase import TBLite

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from ase.optimize.optimize import Optimizer

    from quacc.types import OptParams, OptSchema, RunSchema, VibKwargs, VibThermoSchema


@job
@requires(has_tblite, "tblite must be installed. Refer to the quacc documentation.")
def static_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        xTB method to use
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the TBLite calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = {"method": method}
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc = TBLite(**calc_flags)
    return Recipe(calc).static(atoms, additional_fields=additional_fields)


@job
@requires(has_tblite, "tblite must be installed. Refer to the quacc documentation.")
def relax_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    relax_cell: bool = False,
    fmax: float | None = 0.01,
    max_steps: int = 1000,
    optimizer: type[Optimizer] = BFGS,
    optimizer_kwargs: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    method
        xTB method to use
    fmax
        Maximum force change in eV/A
    max_steps
        Maximum number of steps
    optimizer
        ASE optimizer class to use
    optimizer_kwargs
        Dictionary of keyword arguments to pass to the optimizer
    additional_fields
        Metadata to store in the results
    **calc_kwargs
        Custom kwargs for the tblite calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    _opt_params = calc_kwargs.pop("opt_params", None)  # deprecated
    calc_defaults = {"method": method}
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc = TBLite(**calc_flags)
    return Recipe(calc).relax(
        atoms,
        relax_cell=relax_cell,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        optimizer_kwargs=optimizer_kwargs,
        additional_fields=additional_fields,
        opt_params=_opt_params,
    )


@job
@requires(has_tblite, "tblite must be installed. Refer to the quacc documentation.")
def freq_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    vib_kwargs: VibKwargs | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VibThermoSchema:
    """
    Run a frequency job and calculate thermochemistry.

    Parameters
    ----------
    atoms
        Atoms object
    method
        xTB method to use
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    vib_kwargs
        Dictionary of kwargs for [quacc.runners.ase.Runner.run_vib][].
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the tblite calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.

    Returns
    -------
    VibThermoSchema
        Dictionary of results
    """
    vib_kwargs = vib_kwargs or {}

    calc_defaults = {"method": method}
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc = TBLite(**calc_flags)

    vib_summary = Recipe(calc).vib(
        atoms,
        is_molecule=True,
        vib_kwargs=vib_kwargs,
        additional_fields=additional_fields,
    )
    thermo_summary = ThermoSummarize(
        vib_summary["atoms"],
        vib_summary["results"]["vib_freqs"],
        energy=energy,
        additional_fields=additional_fields,
    ).ideal_gas(temperature=temperature, pressure=pressure)

    return recursive_dict_merge(vib_summary, thermo_summary)
