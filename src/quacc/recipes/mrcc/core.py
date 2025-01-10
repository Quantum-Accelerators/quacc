"""Core recipes for MRCC."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.mrcc._base import run_and_summarize

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "pbe",
    basis: str = "def2-tzvp",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        The method [e.g., PBE or CCSD(T)] to use, this is the value for the calc keyword.
    basis
        Basis set
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize][]
    """
    default_inputs = {"calc": method, "basis": basis, "symm": "off"}

    return run_and_summarize(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        input_swaps=calc_kwargs,
        additional_fields={"name": "MRCC Static"},
        copy_files=copy_files,
    )
