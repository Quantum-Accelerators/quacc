from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize

has_fairchem = bool(find_spec("fairchem"))
has_fairchem_omat = has_fairchem and bool(find_spec("fairchem.data.omat"))

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, SourceDirectory, VaspSchema


@job
@requires(
    has_fairchem_omat,
    "fairchem-data-omat is not installed. Run `pip install quacc[fairchem]`",
)
def omat_static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Carry out a static calculation with OMat settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from fairchem.data.omat.vasp.sets import OMat24StaticSet

    calc_defaults = MPtoASEConverter(atoms=atoms).convert_dict_set(OMat24StaticSet())

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMat Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )
