"""Core recipes for DFTB+"""
from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import covalent as ct
from ase.calculators.dftb import Dftb

from quacc.schemas.ase import summarize_run
from quacc.schemas.atoms import fetch_atoms
from quacc.util.calc import run_calc
from quacc.util.dicts import get_parameters
from quacc.util.files import check_logfile

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import RunSchema

LOG_FILE = "dftb.out"
GEOM_FILE = "geo_end.gen"


@ct.electron
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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    method
        Method to use.
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    atoms = fetch_atoms(atoms)
    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
    }
    flags = get_parameters(defaults, calc_swaps)

    calc = Dftb(**flags)
    final_atoms = run_calc(atoms, calc, geom_file=GEOM_FILE, copy_files=copy_files)

    if check_logfile(LOG_FILE, "SCC is NOT converged"):
        raise ValueError("SCC is not converged")

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "DFTB+ Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    lattice_opt: bool = False,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a structure relaxation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    method
        Method to use.
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    lattice_opt
        Whether to relax the unit cell shape/volume in addition to
        the positions.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    atoms = fetch_atoms(atoms)
    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
        "Driver_": "GeometryOptimization",
        "Driver_LatticeOpt": "Yes" if lattice_opt else "No",
        "Driver_AppendGeometries": "Yes",
        "Driver_MaxSteps": 2000,
    }
    flags = get_parameters(defaults, calc_swaps)

    calc = Dftb(**flags)
    final_atoms = run_calc(atoms, calc, geom_file=GEOM_FILE, copy_files=copy_files)

    if not check_logfile(LOG_FILE, "Geometry converged"):
        raise ValueError("Geometry did not converge")

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "DFTB+ Relax"},
    )
