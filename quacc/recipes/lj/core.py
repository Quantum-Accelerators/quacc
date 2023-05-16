"""Core recipes for EMT"""
from __future__ import annotations

from copy import deepcopy
from typing import Any

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones

from quacc.schemas.calc import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc
from quacc.util.dicts import merge_dicts

# NOTE: This set of minimal recipes is mainly for demonstration purposes


@ct.electron
def static_job(
    atoms: Atoms,
    lj_kwargs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    lj_kwargs
        Dictionary of custom kwargs for the LJ calculator.

    Returns
    -------
    summary
        Summary of the run.
    """

    lj_kwargs = lj_kwargs or {}
    input_atoms = deepcopy(atoms)

    defaults = {"epsilon": 1.0, "sigma": 1.0}
    flags = merge_dicts(defaults, lj_kwargs)

    atoms.calc = LennardJones(**flags)
    atoms = run_calc(atoms)
    summary = summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "LJ Static"}
    )

    return summary


@ct.electron
def relax_job(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    lj_kwargs: dict[str, Any] | None = None,
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
    lj_kwargs
        Dictionary of custom kwargs for the LJ calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    summary
        Summary of the run.
    """

    lj_kwargs = lj_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    defaults = {"epsilon": 1.0, "sigma": 1.0}
    flags = merge_dicts(defaults, lj_kwargs)

    atoms.calc = LennardJones(**flags)
    traj = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )
    summary = summarize_opt_run(
        traj, atoms.calc.parameters, additional_fields={"name": ":J Relax"}
    )

    return summary
