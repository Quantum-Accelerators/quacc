"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT
from ase.optimize import FIRE

from quacc import job
from quacc.runners.calc import run_ase_calc, run_ase_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms,
    calc_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

        !!! Info "Calculator defaults"

            ```python
            {}
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """
    calc_swaps = calc_swaps or {}

    atoms.calc = EMT(**calc_swaps)
    final_atoms = run_ase_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "EMT Static"},
    )


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    calc_swaps: dict[str, Any] | None = None,
    opt_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

        !!! Info "Calculator defaults"

            ```python
            {}
            ```
    opt_swaps
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.calc.run_ase_opt][].

        !!! Info "Optimizer defaults"

            ```python
            {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in
        [quacc.schemas.ase.summarize_opt_run][]
    """
    calc_swaps = calc_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = EMT(**calc_swaps)

    dyn = run_ase_opt(atoms, relax_cell=relax_cell, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
