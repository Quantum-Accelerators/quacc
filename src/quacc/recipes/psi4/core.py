"""Core recipes for Psi4"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.psi4 import Psi4
from monty.dev import requires

from quacc import job
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_run
from quacc.utils.atoms import get_charge_and_spin
from quacc.utils.calc import run_calc
from quacc.utils.dicts import merge_dicts

try:
    import psi4
except ImportError:
    psi4 = None

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import RunSchema


@job
@requires(psi4, "Psi4 not installed. Try conda install -c psi4 psi4")
def static_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    charge
        Charge of the system. If None, this is determined from
        `quacc.utils.atoms.get_charge_and_spin`
    spin_multiplicity
        Multiplicity of the system. If None, this is determined from
        `quacc.utils.atoms.get_charge_and_spin`
    method
        The level of theory to use.
    basis
        Basis set
    calc_swaps
        Dictionary of custom kwargs for the calculator.

        ???+ Note

             Overrides the following defaults:

            ```python
            {
                "mem": "16GB",
                "num_threads": "max",
                "method": method,
                "basis": basis,
                "charge": charge,
                "multiplicity": multiplicity,
                "reference": "uks" if multiplicity > 1 else "rks",
            }
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    charge, spin_multiplicity = get_charge_and_spin(
        atoms, charge=charge, spin_multiplicity=spin_multiplicity
    )

    defaults = {
        "mem": "16GB",
        "num_threads": "max",
        "method": method,
        "basis": basis,
        "charge": charge,
        "multiplicity": spin_multiplicity,
        "reference": "uks" if spin_multiplicity > 1 else "rks",
    }
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Psi4(**flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "Psi4 Static"},
    )
