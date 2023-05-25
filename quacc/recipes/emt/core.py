"""Core recipes for EMT"""
from __future__ import annotations

import covalent as ct
from ase import Atoms
from ase.calculators.emt import EMT

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.atoms import copy_atoms
from quacc.utils.calc import run_ase_opt, run_calc

# NOTE: This set of minimal recipes is mainly for demonstration purposes


@ct.electron
def static_job(atoms: Atoms, emt_kwargs: dict | None = None) -> dict:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    emt_kwargs
        Dictionary of custom kwargs for the EMT calculator.

    Returns
    -------
    summary
        Summary of the run.
    """

    emt_kwargs = emt_kwargs or {}
    input_atoms = copy_atoms(atoms)

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
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    emt_kwargs: dict | None = None,
    opt_kwargs: dict | None = None,
) -> dict:
    """
    Carry out a geometry optimization.

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

    return summarize_opt_run(
        traj, atoms.calc.parameters, additional_fields={"name": "EMT Relax"}
    )
