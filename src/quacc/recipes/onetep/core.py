"""Core recipes for Onetep."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.onetep._base import base_fn

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms, copy_files: list[str] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a basic SCF calculation with ONETEP.

    Parameters
    ----------
    atoms
        The Atoms object.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.onetep.Onetep` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "keywords": {
            "output_detail": "verbose",
            "do_properties": True,
            "cutoff_energy": "600 eV",
            "task": "SinglePoint",
        }
    }

    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "ONETEP Static"},
        copy_files=copy_files,
    )
