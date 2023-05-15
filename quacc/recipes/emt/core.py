"""Core recipes for EMT"""
from __future__ import annotations

from copy import deepcopy
from typing import Any

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.emt import EMT

from quacc.schemas.calc import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc

# NOTE: This set of minimal recipes is mainly for demonstration purposes


@ct.electron
def static_job(
    atoms: Atoms, emt_kwargs: dict[str, Any] | None = None
) -> dict[str, Any]:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    emt_kwargs
        Dictinoary of custom kwargs for the EMT calculator.

    Returns
    -------
    summary
        Summary of the run.
    """

    emt_kwargs = emt_kwargs or {}
    input_atoms = deepcopy(atoms)

    atoms.calc = EMT(**emt_kwargs)
    atoms = run_calc(atoms)
    summary = summarize_run(atoms, input_atoms=input_atoms)

    return summary


@ct.electron
def relax_job(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    emt_kwargs: dict[str, Any] | None = None,
    opt_kwargs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a geometry optimization.

    Parameters
    ----------
    atoms
        .Atoms object
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        .Optimizer class to use for the relaxation.
    emt_kwargs
        Dictionary of custom kwargs for the EMT calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    summary
        Summary of the run.
    """

    emt_kwargs = emt_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    atoms.calc = EMT(**emt_kwargs)
    traj = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )
    summary = summarize_opt_run(traj, atoms.calc.parameters)

    return summary
