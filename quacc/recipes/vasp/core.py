"""Core recipes for VASP"""
from __future__ import annotations

from typing import Any

import covalent as ct
from ase.atoms import Atoms

from quacc.calculators.vasp import Vasp
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts


@ct.electron
def static_job(
    atoms: Atoms, preset: str | None = None, swaps: dict[str, Any] | None = None
) -> dict[str, Any]:
    """
    Function to carry out a single-point calculation.

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
        dictionary of the run summary.
    """

    swaps = swaps or {}

    defaults = {
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lwave": True,
        "nedos": 5001,
        "nsw": 0,
    }
    flags = merge_dicts(defaults, swaps)

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)
    summary = summarize_run(atoms, additional_fields={"name": "VASP Static"})

    return summary


@ct.electron
def relax_job(
    atoms: Atoms,
    preset: str | None = None,
    relax_volume: bool = True,
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to relax a structure.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    relax_volume
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Dictionary of the run summary.
    """

    swaps = swaps or {}

    defaults = {
        "ediffg": -0.02,
        "isif": 3 if relax_volume else 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": False,
        "nsw": 200,
    }
    flags = merge_dicts(defaults, swaps)

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)
    summary = summarize_run(atoms, additional_fields={"name": "VASP Relax"})

    return summary


@ct.electron
def double_relax_job(
    atoms: Atoms,
    preset: str | None = None,
    relax_volume: bool = True,
    swaps1: dict[str, Any] | None = None,
    swaps2: dict[str, Any] | None = None,
) -> dict[dict[str, Any], dict[str, Any]]:
    """
    Function to double-relax a structure. This is particularly useful for
    a few reasons:
    1. To carry out a cheaper pre-relaxation before the high-quality run.
    2. To carry out a GGA calculation before a meta-GGA or hybrid calculation
    that requies the GGA wavefunction.
    3. To carry out volume relaxations where large changes in volume
    can require a second relaxation to resolve forces.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    relax_volume
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
    swaps1
        Dictionary of custom kwargs for the first relaxation.
    swaps2
        Dictionary of custom kwargs for the second relaxation.

    Returns
    -------
    summary
        Dictionary of the run summary.
    """

    swaps1 = swaps1 or {}
    swaps2 = swaps2 or {}

    defaults = {
        "ediffg": -0.02,
        "isif": 3 if relax_volume else 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": True,
        "nsw": 200,
    }

    # Run first relaxation
    flags = merge_dicts(defaults, swaps1)
    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    kpts1 = atoms.calc.kpts
    atoms = run_calc(atoms)
    summary1 = summarize_run(atoms, additional_fields={"name": "VASP DoubleRelax 1"})

    # Run second relaxation
    flags = merge_dicts(defaults, swaps2)
    calc = Vasp(summary1["atoms"], preset=preset, **flags)
    atoms.calc = calc
    kpts2 = atoms.calc.kpts

    # Use ISTART = 0 if this goes from vasp_gam --> vasp_std
    if kpts1 == [1, 1, 1] and kpts2 != [1, 1, 1]:
        atoms.calc.set(istart=0)

    atoms = run_calc(atoms, copy_files=["WAVECAR"])
    summary2 = summarize_run(atoms, additional_fields={"name": "VASP DoubleRelax 2"})

    return {"relax1": summary1, "relax2": summary2}
