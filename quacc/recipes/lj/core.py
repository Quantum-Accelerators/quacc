"""
Core recipes for Lennard-Jones Potential

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""
from __future__ import annotations

from copy import deepcopy

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize import FIRE

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc


@ct.electron
def static_job(
    atoms: Atoms,
    calc_kwargs: dict | None = None,
) -> dict:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    calc_kwargs
        Dictionary of custom kwargs for the LJ calculator

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    calc_kwargs = calc_kwargs or {}
    input_atoms = deepcopy(atoms)

    atoms.calc = LennardJones(**calc_kwargs)
    atoms = run_calc(atoms)

    return summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "LJ Static"}
    )


@ct.electron
def relax_job(
    atoms: Atoms, calc_kwargs: dict | None = None, opt_swaps: dict | None = None
) -> dict:
    """
    Function to carry out a geometry optimization

    Parameters
    ----------
    atoms
        Atoms object
    calc_kwargs
        Dictionary of custom kwargs for the LJ calculator.
    opt_swaps
        Dictionary of swaps for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_opt_run`
    """

    calc_kwargs = calc_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}

    opt_flags = opt_defaults | opt_swaps

    atoms.calc = LennardJones(**calc_kwargs)
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "LJ Relax"})
