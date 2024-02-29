"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import job
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.emt.EMT][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """
    atoms.calc = EMT(**calc_kwargs)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(final_atoms, atoms, additional_fields={"name": "EMT Static"})


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.emt.EMT][] calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """
    opt_params = opt_params or {}

    atoms.calc = EMT(**calc_kwargs)
    dyn = run_opt(atoms, relax_cell=relax_cell, copy_files=copy_files, **opt_params)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
