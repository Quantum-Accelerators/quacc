"""Core recipes for Psi4"""
from __future__ import annotations

from typing import Any, Dict

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.psi4 import Psi4
from monty.dev import requires

try:
    import psi4
except:
    psi4 = None
from quacc.schemas.calc import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc


@ct.electron
@requires(psi4, "Psi4 be installed. Try conda install -c psi4 psi4")
def StaticJob(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    method: str = "wb97x-v",
    basis: str = "def2-tzvp",
    swaps: Dict[str, Any] | None = None,
):
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    charge
        Charge of the system. If None, this is determined from the sum of
        atoms.get_initial_charges().
    mult
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of atoms.get_initial_magnetic_moments().
    method
        The level of theory to use.
    basis
        Basis set
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Summary of the calculation.
    """

    swaps = swaps or {}
    defaults = {
        "mem": "16GB",
        "num_threads": "max",
        "method": method,
        "basis": basis,
        "charge": charge if charge else round(sum(atoms.get_initial_charges())),
        "multiplicity": mult
        if mult
        else round(1 + sum(atoms.get_initial_magnetic_moments())),
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)

    atoms.calc = Psi4(**flags)
    new_atoms = run_calc(atoms)
    summary = summarize_run(new_atoms, input_atoms=atoms)

    return summary
