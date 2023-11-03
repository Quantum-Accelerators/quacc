"""Core recipes for DFTB+"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.dftb import Dftb

from quacc import SETTINGS, job
from quacc.runners.calc import run_ase_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import check_logfile

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase import Atoms

    from quacc.schemas.ase import RunSchema

LOG_FILE = "dftb.out"
GEOM_FILE = "geo_end.gen"


@job
def static_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    calc_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Method to use.
    kpts
        k-point grid to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator that would override the
        calculator defaults. Set a value to `None` to remove a pre-existing key
        entirely. For a list of available keys, refer to the
        `ase.calculators.dftb.Dftb` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
                "Hamiltonian_MaxSccIterations": 200,
                "Hamiltonian_Method": method if "xtb" in method.lower() else None,
                "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
            }
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_MaxSccIterations": 200,
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
    }

    return _base_job(
        atoms,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "DFTB+ Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    relax_cell: bool = False,
    calc_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a structure relaxation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Method to use.
    kpts
        k-point grid to use.
    relax_cell
        Whether to relax the unit cell shape/volume in addition to the
        positions.
    calc_swaps
        Dictionary of custom kwargs for the calculator that would override the
        calculator defaults. Set a value to `None` to remove a pre-existing key
        entirely. For a list of available keys, refer to the
        `ase.calculators.dftb.Dftb` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
                "Hamiltonian_MaxSccIterations": 200,
                "Hamiltonian_Method": method if "xtb" in method.lower() else None,
                "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
                "Driver_": "GeometryOptimization",
                "Driver_LatticeOpt": "Yes" if relax_cell else "No",
                "Driver_AppendGeometries": "Yes", "Driver_MaxSteps": 2000,
            }
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """

    defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_MaxSccIterations": 200,
        "Hamiltonian_Method": method if "xtb" in method.lower() else None,
        "kpts": kpts or ((1, 1, 1) if atoms.pbc.any() else None),
        "Driver_": "GeometryOptimization",
        "Driver_LatticeOpt": "Yes" if relax_cell else "No",
        "Driver_AppendGeometries": "Yes",
        "Driver_MaxSteps": 2000,
    }

    return _base_job(
        atoms,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "DFTB+ Relax"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms,
    defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base job function for DFTB+ recipes.

    Parameters
    ----------
    atoms
        Atoms object
    defaults
        The default calculator parameters to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator that would override the
        calculator defaults. Set a value to `None` to remove a pre-existing key
        entirely. For a list of available keys, refer to the
        `ase.calculators.dftb.Dftb` calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """

    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Dftb(**flags)
    final_atoms = run_ase_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    if SETTINGS.CHECK_CONVERGENCE:
        if check_logfile(LOG_FILE, "SCC is NOT converged"):
            msg = f"SCC is not converged in {LOG_FILE}"
            raise RuntimeError(msg)
        if flags.get("Driver_") == "GeometryOptimization" and not check_logfile(
            LOG_FILE, "Geometry converged"
        ):
            msg = f"Geometry optimization did not complete in {LOG_FILE}"
            raise RuntimeError(msg)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
