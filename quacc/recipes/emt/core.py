"""
Core recipes for EMT

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

from shutil import deepcopy

import covalent as ct
from ase import Atoms
from ase.calculators.emt import EMT

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc


@ct.electron
def static_job(atoms: Atoms, emt_kwargs: dict = None) -> dict:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    emt_kwargs
        Dictionary of custom kwargs for the EMT calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
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
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    emt_kwargs: dict = None,
    opt_kwargs: dict = None,
) -> dict:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
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
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    emt_kwargs = emt_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    atoms.calc = EMT(**emt_kwargs)
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
