"""
QMOF-compatible recipes.

This set of recipes is meant to be compatible with the QMOF Database workflow.
Reference: https://doi.org/10.1016/j.matt.2021.02.015
"""

from __future__ import annotations

from logging import getLogger
from typing import TYPE_CHECKING

from ase.optimize import BFGSLineSearch

from quacc import change_settings, job
from quacc.recipes.vasp._base import run_and_summarize, run_and_summarize_opt

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import OptSchema, SourceDirectory, VaspSchema

    class QMOFRelaxSchema(VaspSchema):
        """Type hint associated with the QMOF relaxation jobs."""

        prerelax_lowacc: VaspSchema | None
        position_relax_lowacc: VaspSchema
        volume_relax_lowacc: VaspSchema | None
        double_relax: VaspSchema


LOGGER = getLogger(__name__)


@job
def qmof_relax_job(
    atoms: Atoms,
    relax_cell: bool = True,
    run_prerelax: bool = True,
    prev_dir: SourceDirectory | None = None,
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
    relax_cell
        True if a volume relaxation should be performed. False if only the
        positions should be updated.
    run_prerelax
        If True, a pre-relax will be carried out with BFGSLineSearch.
        Recommended if starting from hypothetical structures or materials with
        very high starting forces.
    prev_dir
        Source directory copy the WAVECAR from, if present.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    QMOFRelaxSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    # 1. Pre-relaxation
    if run_prerelax:
        summary1 = _prerelax(atoms, **calc_kwargs)
        atoms = summary1["atoms"]
        prev_dir = summary1["dir_name"]

    # 2. Position relaxation (loose)
    summary2 = _loose_relax_positions(atoms, prev_dir=prev_dir, **calc_kwargs)
    atoms = summary2["atoms"]
    prev_dir = summary2["dir_name"]

    # 3. Optional: Volume relaxation (loose)
    if relax_cell:
        summary3 = _loose_relax_cell(atoms, prev_dir=prev_dir, **calc_kwargs)
        atoms = summary3["atoms"]
        prev_dir = summary3["dir_name"]

    # 4. Double Relaxation
    # This is done for two reasons: a) because it can
    # resolve repadding issues when dV is large; b) because we can use LREAL =
    # Auto for the first relaxation and the default LREAL for the second.
    summary4 = _double_relax(
        atoms, prev_dir=prev_dir, relax_cell=relax_cell, **calc_kwargs
    )
    atoms = summary4[-1]["atoms"]
    prev_dir = summary4[-1]["dir_name"]

    # 5. Static Calculation
    summary5 = _static(atoms, prev_dir=prev_dir, **calc_kwargs)
    summary5["prerelax_lowacc"] = summary1 if run_prerelax else None
    summary5["position_relax_lowacc"] = summary2
    summary5["volume_relax_lowacc"] = summary3 if relax_cell else None
    summary5["double_relax"] = summary4

    return summary5


def _prerelax(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> OptSchema:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        Source directory copy the WAVECAR from, if present.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][].
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
    return run_and_summarize_opt(
        atoms,
        preset="QMOFSet",
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults={"fmax": 5.0, "optimizer": BFGSLineSearch},
        additional_fields={"name": "QMOF Prerelax"},
        copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
    )


def _loose_relax_positions(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        Source directory copy the WAVECAR from, if present.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
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
    return run_and_summarize(
        atoms,
        preset="QMOFSet",
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF Loose Relax"},
        copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
    )


def _loose_relax_cell(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        Source directory copy the WAVECAR from, if present.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
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
    return run_and_summarize(
        atoms,
        preset="QMOFSet",
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF Loose Relax Volume"},
        copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
    )


def _double_relax(
    atoms: Atoms,
    prev_dir: SourceDirectory | None = None,
    relax_cell: bool = True,
    **calc_kwargs,
) -> list[VaspSchema]:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        True if a volume relaxation should be performed.
    prev_dir
        Source directory copy the WAVECAR from, if present.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    list[VaspSchema]
        List of dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][]
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

    # To ensure vasp_gam --> vasp_std issues are auto-fixed
    with change_settings({"VASP_USE_CUSTODIAN": True}):
        summary1 = run_and_summarize(
            atoms,
            preset="QMOFSet",
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            additional_fields={"name": "QMOF DoubleRelax 1"},
            copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
        )

    # Update atoms for Relaxation 2
    atoms = summary1["atoms"]

    # Reset LREAL
    del calc_defaults["lreal"]

    # Run second relaxation
    summary2 = run_and_summarize(
        summary1["atoms"],
        preset="QMOFSet",
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF DoubleRelax 2"},
        copy_files={summary1["dir_name"]: ["WAVECAR*"]},
    )
    return [summary1, summary2]


def _static(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Static calculation using production-quality settings.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        Source directory copy the WAVECAR from, if present.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = {
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lwave": True,
        "nsw": 0,
    }
    return run_and_summarize(
        atoms,
        preset="QMOFSet",
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "QMOF Static"},
        copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
    )
