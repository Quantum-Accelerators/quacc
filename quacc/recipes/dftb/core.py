"""Core recipes for DFTB+"""
from __future__ import annotations

import covalent as ct
from ase import Atoms
from ase.calculators.dftb import Dftb

from quacc.schemas.ase import summarize_run
from quacc.util.atoms import copy_atoms
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts
from quacc.util.files import check_logfile

LOG_FILE = "dftb.out"
GEOM_FILE = "geo_end.gen"


@ct.electron
def static_job(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    swaps: dict | None = None,
) -> dict:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Method to use. Accepts "DFTB", "GFN1-xTB", and "GFN2-xTB".
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """

    swaps = swaps or {}
    input_atoms = copy_atoms(atoms)

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
    }
    flags = merge_dicts(defaults, swaps, remove_none=True, auto_lowercase=False)

    atoms.calc = Dftb(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    if check_logfile(LOG_FILE, "SCC is NOT converged"):
        raise ValueError("SCC is not converged")

    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "DFTB+ Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    lattice_opt: bool = False,
    swaps: dict | None = None,
) -> dict:
    """
    Carry out a structure relaxation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Method to use. Accepts 'DFTB', 'GFN1-xTB', and 'GFN2-xTB'.
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    lattice_opt
        Whether to relax the unit cell shape/volume in addition to
        the positions.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """

    swaps = swaps or {}
    input_atoms = copy_atoms(atoms)

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
        "Driver_": "GeometryOptimization",
        "Driver_LatticeOpt": "Yes" if lattice_opt else "No",
        "Driver_AppendGeometries": "Yes",
        "Driver_MaxSteps": 2000,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True, auto_lowercase=False)

    atoms.calc = Dftb(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    if not check_logfile(LOG_FILE, "Geometry converged"):
        raise ValueError("Geometry did not converge")

    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "DFTB+ Relax"},
    )
