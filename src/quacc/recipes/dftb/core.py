"""Core recipes for DFTB+"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.dftb import Dftb

from quacc import job
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_run
from quacc.utils.calc import run_calc
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import check_logfile

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms

    from quacc.schemas.ase import RunSchema

LOG_FILE = "dftb.out"
GEOM_FILE = "geo_end.gen"


@job
def static_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    method
        Method to use.
    kpts
        k-point grid to use. Defaults to None for molecules and (1, 1, 1) for
        solids.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Overrides the following
        defaults:

        ```python
        {
            "Hamiltonian_": "xTB" if "xtb" in method.lower()
            else "DFTB", "Hamiltonian_Method": method if "xtb" in method.lower()
            else None, "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None)
        }
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in `quacc.schemas.ase.RunSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
    }
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Dftb(**flags)
    final_atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    if check_logfile(LOG_FILE, "SCC is NOT converged"):
        msg = "SCC is not converged"
        raise ValueError(msg)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "DFTB+ Static"},
    )


@job
def relax_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    relax_cell: bool = False,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a structure relaxation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    method
        Method to use.
    kpts
        k-point grid to use. Defaults to None for molecules and (1, 1, 1) for
        solids.
    relax_cell
        Whether to relax the unit cell shape/volume in addition to the
        positions.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Overrides the following
        defaults:

        ```python
        {
            "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
            "Hamiltonian_Method": method if "xtb" in method.lower() else None,
            "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
            "Driver_": "GeometryOptimization",
            "Driver_LatticeOpt": "Yes" if relax_cell else "No",
            "Driver_AppendGeometries": "Yes",
            "Driver_MaxSteps": 2000,
        }
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in `quacc.schemas.ase.RunSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
        "Driver_": "GeometryOptimization",
        "Driver_LatticeOpt": "Yes" if relax_cell else "No",
        "Driver_AppendGeometries": "Yes",
        "Driver_MaxSteps": 2000,
    }
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Dftb(**flags)
    final_atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    if not check_logfile(LOG_FILE, "Geometry converged"):
        msg = "Geometry did not converge"
        raise ValueError(msg)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "DFTB+ Relax"},
    )
