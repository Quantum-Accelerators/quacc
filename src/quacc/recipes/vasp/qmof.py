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
from quacc.recipes.vasp.core import _base_job
from quacc.runners.calc import run_ase_opt
from quacc.schemas.ase import summarize_opt_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.ase import OptSchema
    from quacc.schemas.vasp import VaspSchema

    class QMOFRelaxSchema(VaspSchema):
        prerelax_lowacc: VaspSchema | None
        position_relax_lowacc: VaspSchema
        volume_relax_lowacc: VaspSchema | None
        double_relax: VaspSchema


@job
def qmof_relax_job(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    relax_cell: bool = True,
    run_prerelax: bool = True,
    calc_swaps: dict[str, Any] | None = None,
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
        Preset to use from `quacc.calculators.presets.vasp`. Applies for all jobs.
    relax_cell
        True if a volume relaxation should be performed. False if only the
        positions should be updated.
    run_prerelax
        If True, a pre-relax will be carried out with BFGSLineSearch.
        Recommended if starting from hypothetical structures or materials with
        very high starting forces.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. Applies for all jobs.

    Returns
    -------
    QMOFRelaxSchema
        Dictionary of results
    """

    # 1. Pre-relaxation
    if run_prerelax:
        summary1 = _prerelax(atoms, preset, calc_swaps, fmax=5.0)
        atoms = summary1["atoms"]

    # 2. Position relaxation (loose)
    summary2 = _loose_relax_positions(atoms, preset, calc_swaps)
    atoms = summary2["atoms"]

    # 3. Optional: Volume relaxation (loose)
    if relax_cell:
        summary3 = _loose_relax_cell(atoms, preset, calc_swaps)
        atoms = summary3["atoms"]

    # 4. Double Relaxation This is done for two reasons: a) because it can
    # resolve repadding issues when dV is large; b) because we can use LREAL =
    # Auto for the first relaxation and the default LREAL for the second.
    summary4 = _double_relax(atoms, preset, calc_swaps, relax_cell=relax_cell)
    atoms = summary4[1]["atoms"]

    # 5. Static Calculation
    summary5 = _static(atoms, preset, calc_swaps)
    summary5["prerelax_lowacc"] = summary1 if run_prerelax else None
    summary5["position_relax_lowacc"] = summary2
    summary5["volume_relax_lowacc"] = summary3 if relax_cell else None
    summary5["double_relax"] = summary4

    return summary5


def _prerelax(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    calc_swaps: dict[str, Any] | None = None,
    fmax: float = 5.0,
) -> OptSchema:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.
    fmax
        Maximum force in eV/A.

    Returns
    -------
    OptSchema
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    defaults = {
        "auto_kpts": {"kppa": 100},
        "ediff": 1e-4,
        "encut": None,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nelm": 225,
        "nsw": 0,
    }
    flags = merge_dicts(defaults, calc_swaps, remove_nones=False)
    atoms.calc = Vasp(atoms, preset=preset, **flags)
    dyn = run_ase_opt(atoms, fmax=fmax, optimizer=BFGSLineSearch)

    return summarize_opt_run(dyn, additional_fields={"name": "QMOF Prerelax"})


def _loose_relax_positions(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    calc_swaps: dict[str, Any] | None = None,
) -> VaspSchema:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.vasp_summarize_run
    """

    defaults = {
        "auto_kpts": {"kppa": 100},
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
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "QMOF Loose Relax"},
    )


def _loose_relax_cell(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    calc_swaps: dict[str, Any] | None = None,
) -> VaspSchema:
    """
    Volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.vasp_summarize_run
    """

    defaults = {
        "auto_kpts": {"kppa": 100},
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500,
    }
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "QMOF Loose Relax Volume"},
        copy_files=["WAVECAR"],
    )


def _double_relax(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    calc_swaps: dict[str, Any] | None = None,
    relax_cell: bool = True,
) -> list[VaspSchema]:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.
    relax_cell
        True if a volume relaxation should be performed.

    Returns
    -------
    list[VaspSchema]
        List of dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    # Run first relaxation
    defaults = {
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3 if relax_cell else 2,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500 if relax_cell else 250,
    }
    summary1 = _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "QMOF DoubleRelax 1"},
        copy_files=["WAVECAR"],
    )

    # Update atoms for Relaxation 2
    atoms = summary1["atoms"]

    # Reset LREAL
    del defaults["lreal"]

    # Run second relaxation
    summary2 = _base_job(
        summary1["atoms"],
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "QMOF DoubleRelax 2"},
        copy_files=["WAVECAR"],
    )
    return [summary1, summary2]


def _static(
    atoms: Atoms,
    preset: str | None = "QMOFSet",
    calc_swaps: dict[str, Any] | None = None,
) -> VaspSchema:
    """
    Static calculation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.vasp_summarize_run
    """

    defaults = {
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lwave": True,
        "nsw": 0,
    }
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "QMOF Static"},
        copy_files=["WAVECAR"],
    )
