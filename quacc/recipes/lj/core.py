"""Core recipes for Lennard-Jones Potential"""
from __future__ import annotations

import covalent as ct
from ase import Atoms
from ase.calculators.lj import LennardJones

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.atoms import copy_atoms
from quacc.util.calc import run_ase_opt, run_calc
from quacc.util.dicts import merge_dicts

# NOTE: This set of minimal recipes is mainly for demonstration purposes


@ct.electron
def static_job(
    atoms: Atoms,
    lj_kwargs: dict = None,
) -> dict:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    lj_kwargs
        Dictionary of custom kwargs for the LJ calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """

    lj_kwargs = lj_kwargs or {}
    input_atoms = copy_atoms(atoms)

    defaults = {"epsilon": 1.0, "sigma": 1.0}
    flags = merge_dicts(defaults, lj_kwargs)

    atoms.calc = LennardJones(**flags)
    atoms = run_calc(atoms)

    return summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "LJ Static"}
    )


@ct.electron
def relax_job(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    lj_kwargs: dict = None,
    opt_kwargs: dict = None,
) -> dict:
    """
    Function to carry out a geometry optimization.

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
    lj_kwargs
        Dictionary of custom kwargs for the LJ calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    lj_kwargs = lj_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    defaults = {"epsilon": 1.0, "sigma": 1.0}
    flags = merge_dicts(defaults, lj_kwargs)

    atoms.calc = LennardJones(**flags)
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )

    return summarize_opt_run(dyn, additional_fields={"name": "LJ Relax"})
