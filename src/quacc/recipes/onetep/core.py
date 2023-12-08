"""Core recipes for Onetep."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms

from quacc import job
from quacc.calculators.onetep.onetep import OnetepTemplate
from quacc.recipes.onetep._base import base_fn

if TYPE_CHECKING:
    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = None,
    copy_files: list[str] | None = None,
    parallel_info: dict[str] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **calc_kwargs
        Custom kwargs for the Onetep calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.Onetep.Onetep` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
    'pseudo_path': './',
    'output_detail': 'verbose',
    'xc' : 'PBE',
    'do_properties': True,
    'cutoff_energy' : '600 eV',
    'task' : 'SinglePoint'
}

    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "onetep static"},
        copy_files=copy_files,
    )
