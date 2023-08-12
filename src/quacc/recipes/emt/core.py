"""
Core recipes for EMT

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

import warnings

import covalent as ct
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE

from quacc.schemas.ase import DynSchema, RunSchema, summarize_dyn_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc


@ct.electron
def static_job(
    atoms: Atoms | dict,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}

    atoms.calc = EMT(**calc_swaps)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "EMT Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms | dict,
    relax_cell: bool = True,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> DynSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    relax_cell
        Whether to relax the cell
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator
    opt_swaps
        Dictionary of swaps for `run_ase_opt`
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    DynSchema
        Dictionary of results from quacc.schemas.ase.summarize_dyn_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = opt_defaults | opt_swaps

    if relax_cell and not atoms.pbc.any():
        warnings.warn(
            "Volume relaxation requested but no PBCs found. Ignoring.", UserWarning
        )
        relax_cell = False

    atoms.calc = EMT(**calc_swaps)

    if relax_cell:
        atoms = ExpCellFilter(atoms)

    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_dyn_run(dyn, additional_fields={"name": "EMT Relax"})
