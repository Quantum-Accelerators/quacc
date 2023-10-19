"""
Core recipes for EMT

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT
from ase.optimize import FIRE

from quacc import job
from quacc.runners.calc import run_ase_opt, run_calc
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a static calculation.

    ??? Note

        Calculator Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """
    calc_swaps = calc_swaps or {}

    atoms.calc = EMT(**calc_swaps)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "EMT Static"},
    )


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Carry out a geometry optimization.

    ??? Note

        Calculator Defaults:

        ```python
        {}
        ```

        Optimizer Defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator. Overrides the
        following defaults: `{}`
    opt_swaps
        Dictionary of swaps for [quacc.runners.calc.run_ase_opt][].
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][]
    """
    calc_swaps = calc_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = EMT(**calc_swaps)

    dyn = run_ase_opt(atoms, relax_cell=relax_cell, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
