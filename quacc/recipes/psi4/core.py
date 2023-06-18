"""Core recipes for Psi4"""
from __future__ import annotations

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.psi4 import Psi4
from monty.dev import requires

from quacc.schemas.ase import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import remove_dict_empties

try:
    import psi4
except ImportError:
    psi4 = None


@ct.electron
@requires(psi4, "Psi4 not installed. Try conda install -c psi4 psi4")
def static_job(
    atoms: Atoms,
    charge: int | None = None,
    multiplicity: int | None = None,
    method: str = "wb97x-v",
    basis: str = "def2-tzvp",
    swaps: dict | None = None,
) -> dict:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system. If None, this is determined from `atoms.charge`
    multiplicity
        Multiplicity of the system. If None, this is determined from `atoms.spin_multiplicity`
    method
        The level of theory to use.
    basis
        Basis set
    swaps
        Dictionary of custom kwargs for the calculator.
            defaults = {
                "mem": "16GB",
                "num_threads": "max",
                "method": method,
                "basis": basis,
                "charge": charge or getattr(atoms, "charge"),
                "multiplicity": multiplicity or getattr(atoms, "spin_multiplicity"),
                "reference": "uhf" if mult > 1 else None,
            }

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    swaps = swaps or {}

    atoms.charge = charge or getattr(atoms, "charge")
    atoms.spin_multiplicity = multiplicity or getattr(atoms, "spin_multiplicity")

    defaults = {
        "mem": "16GB",
        "num_threads": "max",
        "method": method,
        "basis": basis,
        "charge": atoms.charge,
        "multiplicity": atoms.spin_multiplicity,
        "reference": "uhf" if atoms.spin_multiplicity > 1 else None,
    }
    flags = remove_dict_empties(defaults | swaps)

    atoms.calc = Psi4(**flags)
    new_atoms = run_calc(atoms)

    return summarize_run(
        new_atoms, input_atoms=atoms, additional_fields={"name": "Psi4 Static"}
    )
