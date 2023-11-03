"""Core recipes for VASP."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.calculators.vasp import Vasp
from quacc.runners.calc import run_ase_calc
from quacc.schemas.vasp import vasp_summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.vasp import VaspSchema

    class DoubleRelaxSchema(VaspSchema):
        relax1: VaspSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    calc_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `quacc.calculators.vasp.vasp.Vasp` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "ismear": -5,
                "laechg": True,
                "lcharg": True,
                "lreal": False,
                "lwave": True,
                "nedos": 5001,
                "nsw": 0,
            }
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lwave": True,
        "nedos": 5001,
        "nsw": 0,
    }
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "VASP Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    relax_cell: bool = True,
    calc_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    relax_cell
        True if a volume relaxation (ISIF = 3) should be performed. False if
        only the positions (ISIF = 2) should be updated.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `quacc.calculators.vasp.vasp.Vasp` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "ediffg": -0.02,
                "isif": 3 if relax_cell else 2,
                "ibrion": 2,
                "isym": 0,
                "lcharg": False,
                "lwave": False,
                "nsw": 200,
                "symprec": 1e-8,
            }
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {
        "ediffg": -0.02,
        "isif": 3 if relax_cell else 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": False,
        "nsw": 200,
        "symprec": 1e-8,
    }
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "VASP Relax"},
        copy_files=copy_files,
    )


@job
def double_relax_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    relax_cell: bool = True,
    calc_swaps1: dict[str, Any] | None = None,
    calc_swaps2: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> DoubleRelaxSchema:
    """
    Double-relax a structure. This is particularly useful for a few reasons:

    1. To carry out a cheaper pre-relaxation before the high-quality run.

    2. To carry out a GGA calculation before a meta-GGA or hybrid calculation
    that requires the GGA wavefunction.

    3. To carry out volume relaxations where large changes in volume
    can require a second relaxation to resolve forces.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    relax_cell
        True if a volume relaxation (ISIF = 3) should be performed. False if
        only the positions (ISIF = 2) should be updated.
    calc_swaps1
        Dictionary of custom kwargs for the first relaxation.
    calc_swaps2
        Dictionary of custom kwargs for the second relaxation.
    copy_files
        Files to copy to the (first) runtime directory.

    Returns
    -------
    DoubleRelaxSchema
        Dictionary of results
    """

    # Run first relaxation
    summary1 = relax_job.__wrapped__(
        atoms,
        preset=preset,
        relax_cell=relax_cell,
        calc_swaps=calc_swaps1,
        copy_files=copy_files,
    )

    # Run second relaxation
    summary2 = relax_job.__wrapped__(
        summary1["atoms"],
        preset=preset,
        relax_cell=relax_cell,
        calc_swaps=calc_swaps2,
        copy_files=["WAVECAR"],
    )
    summary2["relax1"] = summary1

    return summary2


def _base_job(
    atoms: Atoms,
    preset: str | None = None,
    defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Base job function for VASP recipes.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.presets.vasp`.
    defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `quacc.calculators.vasp.vasp.Vasp` calculator.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """
    flags = merge_dicts(defaults, calc_swaps, remove_nones=False)

    atoms.calc = Vasp(atoms, preset=preset, **flags)
    atoms = run_ase_calc(atoms, copy_files=copy_files)

    return vasp_summarize_run(atoms, additional_fields=additional_fields)
