"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import job
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory


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
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc = EMT(**calc_kwargs)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return Summarize(additional_fields={"name": "EMT Static"}).run(final_atoms, atoms)


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
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
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
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
    opt_params = opt_params or {}

    calc = EMT(**calc_kwargs)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(
        relax_cell=relax_cell, **opt_params
    )

    return Summarize(additional_fields={"name": "EMT Relax"}).opt(dyn)
