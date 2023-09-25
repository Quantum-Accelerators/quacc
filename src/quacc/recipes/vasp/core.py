"""Core recipes for VASP"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.calculators.vasp import Vasp
from quacc.runners.calc import run_calc
from quacc.schemas import fetch_atoms
from quacc.schemas.vasp import vasp_summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.vasp import VaspSchema

    class DoubleRelaxSchema(VaspSchema):
        relax1: VaspSchema


@job
def static_job(
    atoms: Atoms | dict,
    preset: str | None = None,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Carry out a single-point calculation.

    ??? Note

        Calculator Defaults:

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

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
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
    atoms: Atoms | dict,
    preset: str | None = None,
    relax_cell: bool = True,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Relax a structure.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "ediffg": -0.02,
            "isif": 3 if relax_cell else 2,
            "ibrion": 2,
            "isym": 0,
            "lcharg": False,
            "lwave": False,
            "nsw": 200,
        }
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    relax_cell
        True if a volume relaxation (ISIF = 3) should be performed. False if
        only the positions (ISIF = 2) should be updated.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
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
    atoms: Atoms | dict,
    preset: str | None = None,
    relax_cell: bool = True,
    calc_swaps1: dict | None = None,
    calc_swaps2: dict | None = None,
    copy_files: list[str] | None = None,
) -> DoubleRelaxSchema:
    """
    double_relax a structure. This is particularly useful for a few reasons:

    1. To carry out a cheaper pre-relaxation before the high-quality run.

    2. To carry out a GGA calculation before a meta-GGA or hybrid calculation
    that requires the GGA wavefunction.

    3. To carry out volume relaxations where large changes in volume
    can require a second relaxation to resolve forces.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
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
        summary1,
        preset=preset,
        relax_cell=relax_cell,
        calc_swaps=calc_swaps2,
        copy_files=["WAVECAR"],
    )
    summary2["relax1"] = summary1

    return summary2


def _base_job(
    atoms: Atoms | dict,
    preset: str | None = None,
    defaults: dict | None = None,
    calc_swaps: dict | None = None,
    additional_fields: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Base job function for VASP recipes.


    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """
    atoms = fetch_atoms(atoms)
    flags = merge_dicts(defaults, calc_swaps, remove_empties=False)

    atoms.calc = Vasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms, copy_files=copy_files)

    return vasp_summarize_run(atoms, additional_fields=additional_fields)
