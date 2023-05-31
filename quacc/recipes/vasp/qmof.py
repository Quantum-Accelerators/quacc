"""
QMOF-compatible recipes

his set of recipes is meant to be compatible with the QMOF Database workflow.
Reference: https://doi.org/10.1016/j.matt.2021.02.015
"""
from __future__ import annotations

import covalent as ct
from ase import Atoms

from quacc.calculators.vasp import Vasp
from quacc.schemas.ase import summarize_opt_run
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_ase_opt, run_calc
from quacc.util.dicts import merge_dicts


@ct.electron
def qmof_relax_job(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    relax_volume: bool = True,
    run_prerelax: bool = True,
    swaps: dict = None,
) -> dict:
    """
    Relax a structure in a multi-step process for increased
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
        Atoms object
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
    dict
        Dictionary of results
    """

    swaps = swaps or {}

    # 1. Pre-relaxation
    if run_prerelax:
        summary1 = _prerelax(atoms, preset, swaps, fmax=5.0)
        atoms = summary1["atoms"]

    # 2. Position relaxation (loose)
    summary2 = _loose_relax_positions(atoms, preset, swaps)
    atoms = summary2["atoms"]

    # 3. Optional: Volume relaxation (loose)
    if relax_volume:
        summary3 = _loose_relax_volume(atoms, preset, swaps)
        atoms = summary3["atoms"]

    # 4. Double Relaxation
    # This is done for two reasons: a) because it can resolve repadding
    # issues when dV is large; b) because we can use LREAL = Auto for the
    # first relaxation and the default LREAL for the second.
    summary4 = _double_relax(atoms, preset, swaps, relax_volume=relax_volume)
    atoms = summary4[1]["atoms"]

    # 5. Static Calculation
    summary5 = _static(atoms, preset, swaps)

    return {
        "prerelax-lowacc": summary1 if run_prerelax else None,
        "position-relax-lowacc": summary2,
        "volume-relax-lowacc": summary3 if relax_volume else None,
        "double-relax": summary4,
        "static": summary5,
    }


def _prerelax(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict = None,
    fmax: float = 5.0,
) -> dict:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    fmax
        Maximum force in eV/A.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
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
    dyn = run_ase_opt(atoms, fmax=fmax, optimizer="BFGSLineSearch")

    return summarize_opt_run(dyn, additional_fields={"name": "QMOF Prerelax"})


def _loose_relax_positions(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict = None,
) -> dict:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
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

    return summarize_run(
        atoms, bader=False, additional_fields={"name": "QMOF Loose Relax"}
    )


def _loose_relax_volume(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict = None,
) -> dict:
    """
    Volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
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

    return summarize_run(
        atoms,
        bader=False,
        additional_fields={"name": "QMOF Loose Relax Volume"},
    )


def _double_relax(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict = None,
    relax_volume: bool = True,
) -> dict:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    relax_volume
        True if a volume relaxation should be performed.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
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
    summary1 = summarize_run(
        atoms, bader=False, additional_fields={"name": "QMOF DoubleRelax 1"}
    )
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
    summary2 = summarize_run(
        atoms, bader=False, additional_fields={"name": "QMOF DoubleRelax 2"}
    )

    return [summary1, summary2]


def _static(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    swaps: dict = None,
) -> tuple[Atoms, dict]:
    """
    Static calculation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
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

    return summarize_run(atoms, additional_fields={"name": "QMOF Static"})
