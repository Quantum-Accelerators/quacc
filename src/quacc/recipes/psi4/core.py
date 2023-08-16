"""Core recipes for Psi4"""
from __future__ import annotations

from typing import TYPE_CHECKING

import covalent as ct
from ase.calculators.psi4 import Psi4
from monty.dev import requires

from quacc.schemas.ase import summarize_run
from quacc.schemas.atoms import fetch_atoms
from quacc.util.calc import run_calc
from quacc.util.dicts import remove_dict_empties

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import RunSchema

try:
    import psi4
except ImportError:
    psi4 = None


@ct.electron
@requires(psi4, "Psi4 not installed. Try conda install -c psi4 psi4")
def static_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    multiplicity: int | None = None,
    method: str = "wb97x-v",
    basis: str = "def2-tzvp",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    charge
        Charge of the system. If None, this is determined from the sum of
        `atoms.get_initial_charges()`.
    multiplicity
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of `atoms.get_initial_magnetic_moments()`.
    method
        The level of theory to use.
    basis
        Basis set
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    charge = int(atoms.get_initial_charges().sum()) if charge is None else charge
    multiplicity = (
        int(1 + atoms.get_initial_magnetic_moments().sum())
        if multiplicity is None
        else multiplicity
    )

    defaults = {
        "mem": "16GB",
        "num_threads": "max",
        "method": method,
        "basis": basis,
        "charge": charge,
        "multiplicity": multiplicity,
        "reference": "uks" if multiplicity > 1 else "rks",
    }
    flags = remove_dict_empties(defaults | calc_swaps)

    atoms.calc = Psi4(**flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        charge_and_multiplicity=(charge, multiplicity),
        additional_fields={"name": "Psi4 Static"},
    )
