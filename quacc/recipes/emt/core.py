"""
Core recipes for EMT

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

from copy import deepcopy

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.emt import EMT
from ase.constraints import ExpCellFilter

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc


@ct.electron
def static_job(atoms: Atoms, emt_kwargs: dict | None = None) -> dict:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    emt_kwargs
        Dictionary of custom kwargs for the EMT calculator

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    emt_kwargs = emt_kwargs or {}
    input_atoms = deepcopy(atoms)

    atoms.calc = EMT(**emt_kwargs)
    atoms = run_calc(atoms)

    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "EMT Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms,
    relax_cell: bool = True,
    emt_kwargs: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    emt_kwargs
        Dictionary of custom kwargs for the EMT calculator
    opt_swaps
        Dictionary of swaps for `run_ase_opt`
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "FIRE"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    emt_kwargs = emt_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "FIRE"}
    opt_flags = opt_defaults | opt_swaps

    atoms.calc = EMT(**emt_kwargs)

    if relax_cell:
        atoms = ExpCellFilter(atoms)

    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
