"""Base job functions for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


def base_job(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> RunSchema:
    """
    Base job function for EMT static calculations.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator that would override the
        calculator defaults. Set a value to `None` to remove a pre-existing key
        entirely. For a list of available keys, refer to the
        `ase.calculators.emt.EMT` calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """

    calc_flags = merge_dicts(calc_defaults, calc_swaps)
    atoms.calc = EMT(**calc_flags)
    final_atoms = run_calc(atoms)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields=additional_fields
    )


def base_opt_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_defaults: dict[str, Any] | None = None,
    opt_swaps: dict[str, Any] | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> OptSchema:
    """
    Base job function for EMT geometry optimizations.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    opt_defaults
        The default optimization parameters to use.
    opt_swaps
        Dictionary of custom kwargs for the optimization process that would
        override the optimization defaults. Set a value to `None` to remove a
        pre-existing key entirely. For a list of available keys, refer to
        [quacc.runners.ase.run_opt][].
    calc_defaults
        The default calculator parameters to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator that would override the
        calculator defaults. Set a value to `None` to remove a pre-existing key
        entirely. For a list of available keys, refer to the
        `ase.calculators.dftb.Dftb` calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in
        [quacc.schemas.ase.summarize_opt_run][]
    """

    opt_flags = merge_dicts(opt_defaults, opt_swaps)
    calc_flags = merge_dicts(calc_defaults, calc_swaps)
    atoms.calc = EMT(**calc_flags)

    dyn = run_opt(atoms, relax_cell=relax_cell, **opt_flags)

    return summarize_opt_run(
        dyn, input_atoms=atoms, additional_fields=additional_fields
    )
