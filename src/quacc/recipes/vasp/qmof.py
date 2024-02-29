"""
QMOF-compatible recipes.

This set of recipes is meant to be compatible with the QMOF Database workflow.
Reference: https://doi.org/10.1016/j.matt.2021.02.015
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import BFGSLineSearch

from quacc import job
from quacc.calculators.vasp import Vasp
from quacc.recipes.vasp._base import base_fn
from quacc.runners.ase import run_opt
from quacc.schemas.ase import summarize_opt_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema
    from quacc.schemas._aliases.vasp import QMOFRelaxSchema, VaspSchema


@job
def qmof_relax_job(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    relax_cell: bool = True,
    run_prerelax: bool = True,
    **calc_kwargs,
) -> QMOFRelaxSchema:
    """
    Relax a structure in a multi-step process for increased computational efficiency.
    This is all done in a single compute job. Settings are such that they are compatible
    with the QMOF Database.

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
        Preset to use from `quacc.calculators.vasp.presets`. Applies for all jobs.
    relax_cell
        True if a volume relaxation should be performed. False if only the
        positions should be updated.
    run_prerelax
        If True, a pre-relax will be carried out with BFGSLineSearch.
        Recommended if starting from hypothetical structures or materials with
        very high starting forces.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. Applies for all jobs.

    Returns
    -------
    QMOFRelaxSchema
        Dictionary of results. See the type-hint for the data structure.
    """

    # 1. Pre-relaxation
    if run_prerelax:
        summary1 = _prerelax(atoms, preset, fmax=5.0, **calc_kwargs)
        atoms = summary1["atoms"]

    # 2. Position relaxation (loose)
    summary2 = _loose_relax_positions(atoms, preset, **calc_kwargs)
    atoms = summary2["atoms"]

    # 3. Optional: Volume relaxation (loose)
    if relax_cell:
        summary3 = _loose_relax_cell(atoms, preset, **calc_kwargs)
        atoms = summary3["atoms"]

    # 4. Double Relaxation This is done for two reasons: a) because it can
    # resolve repadding issues when dV is large; b) because we can use LREAL =
    # Auto for the first relaxation and the default LREAL for the second.
    summary4 = _double_relax(atoms, preset, relax_cell=relax_cell, **calc_kwargs)
    atoms = summary4[1]["atoms"]

    # 5. Static Calculation
    summary5 = _static(atoms, preset, **calc_kwargs)
    summary5["prerelax_lowacc"] = summary1 if run_prerelax else None
    summary5["position_relax_lowacc"] = summary2
    summary5["volume_relax_lowacc"] = summary3 if relax_cell else None
    summary5["double_relax"] = summary4

    return summary5


def _prerelax(
    atoms: Atoms, preset: str | None = "QMOFSet", fmax: float = 5.0, **calc_kwargs
) -> OptSchema:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    fmax
        Maximum force in eV/A.
    **kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.
    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "pmg_kpts": {"kppa": 100},
        "ediff": 1e-4,
        "encut": None,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nelm": 225,
        "nsw": 0,
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    atoms.calc = Vasp(atoms, preset=preset, **calc_flags)
    dyn = run_opt(atoms, fmax=fmax, optimizer=BFGSLineSearch)

    return summarize_opt_run(dyn, additional_fields={"name": "QMOF Prerelax"})


def _loose_relax_positions(
    atoms: Atoms, preset: str | None = "QMOFSet", **calc_kwargs
) -> VaspSchema:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    **kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "pmg_kpts": {"kppa": 100},
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
    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF Loose Relax"},
    )


def _loose_relax_cell(
    atoms: Atoms, preset: str | None = "QMOFSet", **calc_kwargs
) -> VaspSchema:
    """
    Volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "pmg_kpts": {"kppa": 100},
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500,
    }
    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF Loose Relax Volume"},
    )


def _double_relax(
    atoms: Atoms, preset: str | None = "QMOFSet", relax_cell: bool = True, **calc_kwargs
) -> list[VaspSchema]:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    relax_cell
        True if a volume relaxation should be performed.
    **calc_kwargs
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.
    Returns
    -------
    list[VaspSchema]
        List of dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
        See the type-hint for the data structure.
    """

    # Run first relaxation
    calc_defaults = {
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3 if relax_cell else 2,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500 if relax_cell else 250,
    }
    summary1 = base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF DoubleRelax 1"},
    )

    # Update atoms for Relaxation 2
    atoms = summary1["atoms"]

    # Reset LREAL
    del calc_defaults["lreal"]

    # Run second relaxation
    summary2 = base_fn(
        summary1["atoms"],
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF DoubleRelax 2"},
    )
    return [summary1, summary2]


def _static(atoms: Atoms, preset: str | None = "QMOFSet", **calc_kwargs) -> VaspSchema:
    """
    Static calculation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    **kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lwave": True,
        "nsw": 0,
    }
    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF Static"},
    )
