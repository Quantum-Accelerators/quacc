"""Core recipes for EMT"""
from __future__ import annotations

from copy import deepcopy
from typing import Any

from ase.atoms import Atoms
from ase.calculators.emt import EMT
from ase.io import read
from ase.optimize import FIRE

from quacc.schemas.calc import summarize_opt_run, summarize_run

# NOTE: This set of minimal recipes is mainly for demonstration purposes


def static_job(atoms: Atoms, asap_cutoff: bool = False) -> dict[str, Any]:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    asap_cutoff
        If an ASAP-style cutoff should be used.

    Returns
    -------
    summary
        Summary of the run.
    """
    input_atoms = deepcopy(atoms)

    atoms.calc = EMT(asap_cutoff=asap_cutoff)
    atoms.get_potential_energy()
    summary = summarize_run(atoms, input_atoms=input_atoms)

    return summary


def relax_job(
    atoms: Atoms,
    asap_cutoff: bool = False,
    fmax: float = 0.03,
    opt_kwargs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a geometry optimization.

    Parameters
    ----------
    atoms
        .Atoms object
    asap_cutoff
        If an ASAP-style cutoff should be used.
    fmax
        Tolerance for the force convergence (in eV/A).
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    summary
        Summary of the run.
    """
    if opt_kwargs is None:
        opt_kwargs = {}

    atoms.calc = EMT(asap_cutoff=asap_cutoff)
    dyn = FIRE(atoms, logfile="opt.log", trajectory="opt.traj", **opt_kwargs)
    dyn.run(fmax=fmax)
    traj = read("opt.traj", index=":")
    summary = summarize_opt_run(traj, atoms.calc.parameters)

    return summary
