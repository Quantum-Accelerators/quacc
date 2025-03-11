"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.calculators.calculator import BaseCalculator

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    calc: BaseCalculator,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    calc
        Instantiated ASE calculator.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()
    return Summarize(additional_fields=additional_fields).run(final_atoms, atoms)


@job
def relax_job(
    atoms: Atoms,
    calc: BaseCalculator,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
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

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    opt_params = opt_params or {}
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(
        relax_cell=relax_cell, **opt_params
    )

    return Summarize(additional_fields=additional_fields).opt(dyn)
