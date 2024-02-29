"""Core recipes for Psi4."""

from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.recipes.psi4._base import base_fn

try:
    import psi4
except ImportError:
    psi4 = None

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
@requires(psi4, "Psi4 not installed. Try conda install -c psi4 psi4")
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97x-v",
    basis: str = "def2-tzvp",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        The level of theory to use.
    basis
        Basis set
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **kwargs
        Custom kwargs for the Psi4 calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.psi4.Psi4][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "mem": "16GB",
        "num_threads": "max",
        "method": method,
        "basis": basis,
        "charge": charge,
        "multiplicity": spin_multiplicity,
        "reference": "uks" if spin_multiplicity > 1 else "rks",
    }
    return base_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=kwargs,
        additional_fields={"name": "Psi4 Static"},
        copy_files=copy_files,
    )
