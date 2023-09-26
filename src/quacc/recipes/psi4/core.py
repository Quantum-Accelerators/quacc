"""Core recipes for Psi4"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.psi4 import Psi4
from monty.dev import requires

from quacc import job
from quacc.runners.calc import run_calc
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_run
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
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97x-v",
    basis: str = "def2-tzvp",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "mem": "16GB",
            "num_threads": "max",
            "method": method,
            "basis": basis,
            "charge": charge,
            "multiplicity": spin_multiplicity,
            "reference": "uks" if spin_multiplicity > 1 else "rks",
        }
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        The level of theory to use.
    basis
        Basis set
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    defaults = {
        "mem": "16GB",
        "num_threads": "max",
        "method": method,
        "basis": basis,
        "charge": charge,
        "multiplicity": spin_multiplicity,
        "reference": "uks" if spin_multiplicity > 1 else "rks",
    }
    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "Psi4 Static"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms | dict,
    charge: int,
    spin_multiplicity: int,
    defaults: dict | None = None,
    calc_swaps: dict | None = None,
    additional_fields: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base function to carry out Psi4 recipes.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    defaults
        The default calculator parameters.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms = fetch_atoms(atoms)
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Psi4(**flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )
