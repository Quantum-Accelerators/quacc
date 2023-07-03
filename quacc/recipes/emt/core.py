"""
Core recipes for EMT

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

import warnings
from typing import Literal

import covalent as ct
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc


@ct.electron
def static_job(
    input_atoms: Atoms | dict[Literal["atoms"], Atoms], calc_kwargs: dict | None = None
) -> dict:
    """
    Carry out a static calculation.

    Parameters
    ----------
    input_atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    calc_kwargs
        Dictionary of custom kwargs for the EMT calculator

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    input_atoms = input_atoms["atoms"] if isinstance(input_atoms, dict) else input_atoms
    calc_kwargs = calc_kwargs or {}

    input_atoms.calc = EMT(**calc_kwargs)
    atoms = run_calc(input_atoms)

    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "EMT Static"},
    )


@ct.electron
def relax_job(
    input_atoms: Atoms | dict[Literal["atoms"], Atoms],
    relax_cell: bool = True,
    calc_kwargs: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    input_atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    relax_cell
        Whether to relax the cell
    calc_kwargs
        Dictionary of custom kwargs for the EMT calculator
    opt_swaps
        Dictionary of swaps for `run_ase_opt`
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """
    input_atoms = input_atoms["atoms"] if isinstance(input_atoms, dict) else input_atoms
    calc_kwargs = calc_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = opt_defaults | opt_swaps

    if relax_cell and not input_atoms.pbc.any():
        warnings.warn(
            "Volume relaxation requested but no PBCs found. Ignoring.", UserWarning
        )
        relax_cell = False

    input_atoms.calc = EMT(**calc_kwargs)

    if relax_cell:
        input_atoms = ExpCellFilter(input_atoms)

    dyn = run_ase_opt(input_atoms, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
