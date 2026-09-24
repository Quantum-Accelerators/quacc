"""Recipes to reproduce OMol settings"""

from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING, Any

from monty.dev import requires

from quacc import job

if TYPE_CHECKING:
    from typing import Literal

    from ase.atoms import Atoms
    from fairchem.data.omol.orca.calc import Vertical

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
    vertical: Vertical | None = None,
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | Literal["max"] = 12,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
):
    """
    Carry out a single-point calculation using the OMol settings.

    The custom `def2-tzvpd.bas` basis set file shipped with fairchem is
    automatically copied into the runtime directory.

    Parameters
    ----------
    atoms
        Atoms object.
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    vertical
        Vertical of the OMol dataset to use. Defaults to `Vertical.Default`.
    orcasimpleinput
        List of `orcasimpleinput` keywords. If provided, this **replaces** (rather
        than updates) fairchem's default OMol keywords (e.g. RIJCOSX, def2/J,
        DEFGRID3, NormalConv), so the calculation will no longer be OMol-compatible
        unless those keywords are included. The functional, basis set, and EnGrad
        are always included.
    orcablocks
        List of `orcablocks` entries. If provided, this **replaces** (rather than
        updates) fairchem's default OMol blocks (e.g. SCF convergence thresholds
        and the custom basis set file), so the calculation will no longer be
        OMol-compatible unless those blocks are included.
    nprocs
        Number of processors to use in the `%pal` block. Use "max" to use
        all physical cores.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.

    Returns
    -------
    dict
        Dictionary of results.
    """
    import fairchem.data.omol.orca as omol_orca
    from fairchem.data.omol.orca.calc import Vertical
    from fairchem.data.omol.orca.recipes import single_point_calculation

    basis_file = {
        "source": Path(omol_orca.__file__).parent / "basis",
        "filenames": ["def2-tzvpd.bas"],
    }

    return single_point_calculation(
        atoms,
        charge,
        spin_multiplicity,
        vertical=Vertical.Default if vertical is None else vertical,
        orcasimpleinput=None if orcasimpleinput is None else list(orcasimpleinput),
        orcablocks=None if orcablocks is None else list(orcablocks),
        nprocs=nprocs,
        copy_files=[basis_file, *(copy_files or [])],
        additional_fields=additional_fields,
    )
