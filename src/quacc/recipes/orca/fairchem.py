"""Recipes to reproduce OMol settings"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING, Any

from fairchem.data.omol.orca.calc import Vertical
from fairchem.data.omol.orca.recipes import single_point_calculation
from monty.dev import requires

from quacc import job

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import CopyFiles

has_fairchem = bool(find_spec("fairchem"))
has_fairchem_omol = (
    has_fairchem
    and bool(find_spec("fairchem.data"))
    and bool(find_spec("fairchem.data.omol"))
)


@job
@requires(
    has_fairchem_omol,
    "fairchem-data-omol is not installed. Run `pip install quacc[fairchem]`",
)
def omol_static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    vertical: Vertical = Vertical.Default,
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
):
    """
    Carry out a single-point calculation using the OMol settings.

    Parameters
    ----------
    atoms
        Atoms object.
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    vertical
        Vertical of the OMol dataset to use.
    orcasimpleinput
        List of `orcasimpleinput` swaps for the calculator.
    orcablocks
        List of `orcablocks` swaps for the calculator.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.

    Returns
    -------
    dict
        Dictionary of results.
    """
    return single_point_calculation(
        atoms,
        charge,
        spin_multiplicity,
        vertical=vertical,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
        copy_files=copy_files,
        additional_fields=additional_fields,
    )
