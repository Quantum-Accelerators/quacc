"""
Core recipes for the Q-Chem
"""
from __future__ import annotations

from copy import deepcopy
from typing import Literal

from ase.atoms import Atoms
from monty.dev import requires

from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo_run,
    summarize_vib_run,
)
from quacc.util.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.util.thermo import ideal_gas
from quacc.calculators.qchem import QChem

def static_job(
    atoms: Atoms,
    qchem_kwargs: dict,
) -> dict:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    qchem_kwargs
        Dictionary of kwargs for the Q-Chem calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """
    input_atoms = deepcopy(atoms)

    calc = QChem(atoms, **qchem_kwargs)
    atoms.calc = calc
    atoms = run_calc(atoms)
    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "Q-Chem Static"},
    )


def opt_job(
    atoms: Atoms,
    qchem_kwargs: dict,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "Sella",
    opt_kwargs: dict | None = None,
) -> dict:
    """
    Optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object
    qchem_kwargs
        Dictionary of kwargs for the Q-Chem calculator.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        .Optimizer class to use for the optimization.
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    opt_kwargs = opt_kwargs or {}
    if optimizer == "Sella":
        if "internal" not in opt_kwargs:
            opt_kwargs["internal"] = True
        if "order" not in opt_kwargs:
            opt_kwargs["order"] = 0

    calc = QChem(atoms, **qchem_kwargs)
    atoms.calc = calc
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )

    return summarize_opt_run(dyn, additional_fields={"name": "Q-Chem Optimization"})

