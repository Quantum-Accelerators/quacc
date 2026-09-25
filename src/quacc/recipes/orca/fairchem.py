"""Recipes to reproduce OMol settings"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING, Any

import psutil
from fairchem.data.omol.orca.calc import ORCA_BLOCKS, ORCA_SIMPLE_INPUT, Vertical
from fairchem.data.omol.orca.recipes import single_point_calculation
from monty.dev import requires

from quacc import job
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Literal

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
    nprocs: int | Literal["max"] = "max",
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
        List of `orcasimpleinput` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcasimpleinput`.
    orcablocks
        List of `orcablocks` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcablocks`.
    nprocs
        Number of processors to use. Defaults to the number of physical cores.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.

    Returns
    -------
    dict
        Dictionary of results.
    """
    nprocs = (psutil.cpu_count(logical=False) if nprocs == "max" else nprocs) or 1

    return single_point_calculation(
        atoms,
        charge,
        spin_multiplicity,
        vertical=vertical,
        orcasimpleinput=merge_list_params(ORCA_SIMPLE_INPUT, orcasimpleinput)
        if orcasimpleinput
        else None,
        orcablocks=merge_list_params(ORCA_BLOCKS, orcablocks) if orcablocks else None,
        nprocs=nprocs,
        copy_files=copy_files,
        additional_fields=additional_fields,
    )
