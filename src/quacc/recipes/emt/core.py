"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import job
from quacc.recipes.common.core import relax_job as relax_job_
from quacc.recipes.common.core import static_job as static_job_

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms, additional_fields: dict[str, Any] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.emt.EMT][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    return static_job_(
        atoms,
        EMT(**calc_kwargs),
        additional_fields={"name": "EMT Static"} | (additional_fields or {}),
    )


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.emt.EMT][] calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    return relax_job_(
        atoms,
        EMT(**calc_kwargs),
        relax_cell=relax_cell,
        opt_params=opt_params,
        copy_files=copy_files,
        additional_fields={"name": "EMT Relax"} | (additional_fields or {}),
    )
