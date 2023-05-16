"""Core recipes for DFTB+"""
from __future__ import annotations

from copy import deepcopy
from shutil import which
from typing import Any

import covalent as ct
from ase import Atoms
from ase.calculators.dftb import Dftb
from monty.dev import requires

from quacc.schemas.calc import summarize_run
from quacc.util.calc import _check_logfile, run_calc
from quacc.util.dicts import merge_dicts

DFTBPLUS_EXISTS = bool(which("dftb+"))
LOG_FILE = "dftb.out"
GEOM_FILE = "geo_end.gen"


@requires(
    DFTBPLUS_EXISTS,
    "DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
@ct.electron
def static_job(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict[str, Any] | None = None,
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    method
        Method to use. Accepts 'DFTB', 'GFN1-xTB', and 'GFN2-xTB'.
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    atoms
        .Atoms object after the run.
    summary
        Dictionary of results from the calculation.
    """

    swaps = swaps or {}
    input_atoms = deepcopy(atoms)

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts if kpts else (1, 1, 1) if atoms.pbc.any() else None,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True, auto_lowercase=False)

    atoms.calc = Dftb(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)
    scc_check = _check_logfile(LOG_FILE, "SCC is NOT converged")
    if scc_check:
        raise ValueError("SCC is not converged")
    summary = summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "DFTB+ Static"}
    )

    return summary


@requires(
    DFTBPLUS_EXISTS,
    "DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
@ct.electron
def relax_job(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict[str, Any] | None = None,
    lattice_opt: bool = False,
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a structure relaxation.

    Parameters
    ----------
    atoms
        .Atoms object
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
    summary
        Dictionary of results from the calculation.
    """

    swaps = swaps or {}
    input_atoms = deepcopy(atoms)

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts if kpts else (1, 1, 1) if atoms.pbc.any() else None,
        "Driver_": "GeometryOptimization",
        "Driver_LatticeOpt": "Yes" if lattice_opt else "No",
        "Driver_AppendGeometries": "Yes",
        "Driver_MaxSteps": 2000,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True, auto_lowercase=False)

    atoms.calc = Dftb(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)
    geom_check = _check_logfile(LOG_FILE, "Geometry converged")
    if not geom_check:
        raise ValueError("Geometry did not converge")
    summary = summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "DFTB+ Relax"}
    )

    return summary
