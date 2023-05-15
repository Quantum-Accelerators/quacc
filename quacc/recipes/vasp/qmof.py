"""QMOF-compatible recipes"""
from __future__ import annotations

from typing import Any

import covalent as ct
from ase.atoms import Atoms

from quacc.calculators.vasp import Vasp
from quacc.schemas.calc import summarize_opt_run
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_ase_opt, run_calc
from quacc.util.dicts import merge_dicts

# This set of recipes is meant to be compatible with the QMOF Database workflow.
# Reference: https://doi.org/10.1016/j.matt.2021.02.015


def qmof_relax_job(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    relax_volume: bool = True,
    run_prerelax: bool = True,
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to relax a structure in a multi-step process for increased
    computational efficiency. This is all done in a single compute job.
    Settings are such that they are compatible with the QMOF Database.

    1. A "pre-relaxation" with BFGSLineSearch to resolve very high forces.
    2. Position relaxation with default ENCUT and coarse k-point grid.
    3. Optional: volume relaxation with coarse k-point grid.
    4. Double relaxation using production-quality settings.
    5. Static calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use. Applies for all jobs.
    relax_volume
        True if a volume relaxation should be performed.
        False if only the positions should be updated.
    run_prerelax
        If True, a pre-relax will be carried out with BFGSLineSearch.
        Recommended if starting from hypothetical structures or materials
        with very high starting forces.
    swaps
        Dictionary of custom kwargs for the calculator. Applies for all jobs.

    Returns
    -------
    summary
        Dictionary of the run summary.
    """

    swaps = swaps or {}

    # 1. Pre-relaxation
    if run_prerelax:
        summary1 = prerelax(atoms, preset, swaps, fmax=5.0)
        atoms = summary1["atoms"]

    # 2. Position relaxation (loose)
    summary2 = loose_relax_positions(atoms, preset, swaps)
    atoms = summary2["atoms"]

    # 3. Optional: Volume relaxation (loose)
    if relax_volume:
        summary3 = loose_relax_volume(atoms, preset, swaps)
        atoms = summary3["atoms"]

    # 4. Double Relaxation
    # This is done for two reasons: a) because it can resolve repadding
    # issues when dV is large; b) because we can use LREAL = Auto for the
    # first relaxation and the default LREAL for the second.
    summary4 = double_relax(atoms, preset, swaps, relax_volume=relax_volume)
    atoms = summary4["relax2"]["atoms"]

    # 5. Static Calculation
    summary5 = static(atoms, preset, swaps)

    return {
        "prerelax-lowacc": summary1 if prerelax else None,
        "position-relax-lowacc": summary2,
        "volume-relax-lowacc": summary3 if relax_volume else None,
        "double-relax": summary4,
        "static": summary5,
    }


@ct.electron
def prerelax(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict[str, Any] | None = None,
    fmax: float = 5.0,
) -> dict[str, Any]:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    fmax
        Maximum force in eV/A.

    Returns
    -------
    summary
        Summary of the run.
    """

    swaps = swaps or {}

    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediff": 1e-4,
        "encut": None,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nelm": 225,
        "nsw": 0,
    }
    flags = merge_dicts(defaults, swaps)
    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    traj = run_ase_opt(atoms, fmax=fmax, optimizer="BFGSLineSearch")

    summary = summarize_opt_run(traj, calc.parameters)

    return summary


@ct.electron
def loose_relax_positions(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Summary of the run.
    """

    swaps = swaps or {}

    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediff": 1e-4,
        "ediffg": -0.05,
        "encut": None,
        "ibrion": 2,
        "isif": 2,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 250,
    }
    flags = merge_dicts(defaults, swaps)
    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)

    summary = summarize_run(atoms, bader=False)

    return summary


@ct.electron
def loose_relax_volume(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Summary of the run.
    """

    swaps = swaps or {}

    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500,
    }
    flags = merge_dicts(defaults, swaps)
    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms, copy_files=["WAVECAR"])

    summary = summarize_run(atoms, bader=False)

    return summary


@ct.electron
def double_relax(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict[str, Any] | None = None,
    relax_volume: bool = True,
) -> dict[str, Any]:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    relax_volume
        True if a volume relaxation should be performed.

    Returns
    -------
    summary
        Summary of the run.
    """

    swaps = swaps or {}

    defaults = {
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3 if relax_volume else 2,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500 if relax_volume else 250,
    }

    # Run first relaxation
    flags = merge_dicts(defaults, swaps)
    calc1 = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc1
    atoms = run_calc(atoms, copy_files=["WAVECAR"])

    # Update atoms for
    summary1 = summarize_run(atoms, bader=False)
    atoms = summary1["atoms"]

    # Reset LREAL
    del defaults["lreal"]

    # Run second relaxation
    flags = merge_dicts(defaults, swaps)
    calc2 = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc2

    # Use ISTART = 0 if this goes from vasp_gam --> vasp_std
    if calc1.kpts == [1, 1, 1] and calc2.kpts != [1, 1, 1]:
        atoms.calc.set(istart=0)

    atoms = run_calc(atoms, copy_files=["WAVECAR"])
    summary2 = summarize_run(atoms, bader=False)

    return {"relax1": summary1, "relax2": summary2}


@ct.electron
def static(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict[str, Any] | None = None,
) -> tuple[Atoms, dict[str, Any]]:
    """
    Static calculation using production-quality settings.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Summary of the run.
    """

    swaps = swaps or {}

    defaults = {
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lwave": True,
        "nsw": 0,
    }

    # Run static calculation
    flags = merge_dicts(defaults, swaps)
    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms, copy_files=["WAVECAR"])

    summary = summarize_run(atoms)

    return summary
