"""Core recipes for Q-Chem."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.mrcc._base import run_and_summarize

has_sella = bool(find_spec("sella"))

if has_sella:
    pass

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "pbe",
    basis: str = "sto-3g",
    mrccinput: dict[str, str] | None = None,
    mrccblocks: str | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    mrccinput
        Dictionary of `mrccinput` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.mrcc.MRCC][] calculator for details on `mrccinput`.
    mrccblocks
        String for the `mrccblocks` input.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][].
        See the type-hint for the data structure.
    """
    default_inputs = {"calc": method, "basis": basis}

    return run_and_summarize(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        blocks=mrccblocks,
        input_swaps=mrccinput,
        additional_fields={"name": "MRCC Static"},
        copy_files=copy_files,
    )
