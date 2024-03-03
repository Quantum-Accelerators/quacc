"""Schemas for VASP."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING

from ase.io import read
from emmet.core.tasks import TaskDoc
from monty.os.path import zpath
from pymatgen.command_line.bader_caller import bader_analysis_from_path
from pymatgen.command_line.chargemol_caller import ChargemolAnalysis
from pymatgen.entries.compatibility import (
    CompatibilityError,
    MaterialsProject2020Compatibility,
)

from quacc import SETTINGS
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import clean_task_doc
from quacc.wflow_tools.db import results_to_db

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from maggma.core import Store

    from quacc.schemas._aliases.vasp import BaderSchema, ChargemolSchema, VaspSchema

logger = logging.getLogger(__name__)


def vasp_summarize_run(
    final_atoms: Atoms,
    dir_path: str | Path | None = None,
    move_magmoms: bool = True,
    run_bader: bool | None = None,
    run_chargemol: bool | None = None,
    check_convergence: bool | None = None,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    store: Store | None = None,
) -> VaspSchema:
    """
    Get tabulated results from a VASP run and store them in a database-friendly format.

    Parameters
    ----------
    final_atoms
        ASE Atoms object following a calculation.
    dir_path
        Path to VASP outputs. A value of None specifies the calculator directory.
    move_magmoms
        Whether to move the final magmoms of the original Atoms object to the
        initial magmoms of the returned Atoms object.
    run_bader
        Whether a Bader analysis should be performed. Will not run if bader
        executable is not in PATH even if bader is set to True. Defaults to
        VASP_BADER in settings.
    run_chargemol
        Whether a Chargemol analysis should be performed. Will not run if chargemol
        executable is not in PATH even if chargmeol is set to True. Defaults to
        VASP_CHARGEMOL in settings.
    check_convergence
        Whether to throw an error if convergence is not reached. Defaults to True in
        settings.
    report_mp_corrections
        Whether to apply the MP corrections to the task document. Defaults to False.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.STORE` will be used.

    Returns
    -------
    VaspSchema
        Dictionary representation of the task document
    """

    additional_fields = additional_fields or {}
    run_bader = SETTINGS.VASP_BADER if run_bader is None else run_bader
    run_chargemol = SETTINGS.VASP_CHARGEMOL if run_chargemol is None else run_chargemol
    check_convergence = (
        SETTINGS.CHECK_CONVERGENCE if check_convergence is None else check_convergence
    )
    dir_path = Path(dir_path or final_atoms.calc.directory)
    store = SETTINGS.STORE if store is None else store

    # Fetch all tabulated results from VASP outputs files. Fortunately, emmet
    # already has a handy function for this
    vasp_task_model = TaskDoc.from_directory(dir_path)

    # Get MP corrections
    if report_mp_corrections:
        mp_compat = MaterialsProject2020Compatibility()
        try:
            corrected_entry = mp_compat.process_entry(
                vasp_task_model.structure_entry, on_error="raise"
            )
            vasp_task_model.entry = corrected_entry
        except CompatibilityError as err:
            logger.warning(err)

    # Convert the VASP task model to a dictionary
    vasp_task_doc = vasp_task_model.model_dump()

    # Check for calculation convergence
    if check_convergence and vasp_task_doc["state"] != "successful":
        raise RuntimeError(
            f"VASP calculation did not converge. Will not store task data. Refer to {dir_path}"
        )

    initial_atoms = read(zpath(dir_path / "POSCAR"))
    base_task_doc = summarize_run(
        final_atoms, initial_atoms, move_magmoms=move_magmoms, store=False
    )

    # Get Bader analysis
    if run_bader:
        try:
            bader_results = _bader_runner(dir_path)
        except Exception:
            bader_results = None
            logging.warning("Bader analysis could not be performed.", exc_info=True)

        if bader_results:
            vasp_task_doc["bader"] = bader_results

    # Get the Chargemol analysis
    if run_chargemol:
        try:
            chargemol_results = _chargemol_runner(dir_path)
        except Exception:
            chargemol_results = None
            logging.warning("Chargemol analysis could not be performed.", exc_info=True)

        if chargemol_results:
            vasp_task_doc["chargemol"] = chargemol_results

    # Make task document
    unsorted_task_doc = vasp_task_doc | base_task_doc | additional_fields
    task_doc = clean_task_doc(unsorted_task_doc)

    # Store the results
    if store:
        results_to_db(store, task_doc)

    return task_doc


def _bader_runner(path: Path | str) -> BaderSchema:
    """
    Runs a Bader partial charge and spin moment analysis using the VASP output
    files in the given path. This function requires that `bader` is located in
    your PATH environment variable. See
    http://theory.cm.utexas.edu/henkelman/code/bader for the bader code.

    Parameters
    ----------
    path
        The path where the VASP output files are located. Must include CHGCAR,
        AECCAR0, AECCAR2, and POTCAR files. These files can be gzip'd or not --
        it doesn't matter. If None, the current working directory is used.
    structure
        The structure object to attach the Bader charges and spins to.

    Returns
    -------
    BaderSchema
        Dictionary containing the Bader analysis summary
    """

    # Make sure files are present
    relevant_files = ["AECCAR0", "AECCAR2", "CHGCAR", "POTCAR"]
    for f in relevant_files:
        if not Path(path, f).exists() and not Path(path, f"{f}.gz").exists():
            msg = f"Could not find {f} in {path}."
            raise FileNotFoundError(msg)

    bader_stats = bader_analysis_from_path(path)

    # Store the partial charge, which is much more useful than the raw charge
    # and is more intuitive than the charge transferred. An atom with a positive
    # partial charge is cationic, whereas an atom with a negative partial charge
    # is anionic.
    bader_stats["partial_charges"] = [-c for c in bader_stats["charge_transfer"]]

    # Some cleanup of the returned dictionary
    if "magmom" in bader_stats:
        bader_stats["spin_moments"] = bader_stats["magmom"]
    for k in ["charge", "charge_transfer", "reference_used", "magmom"]:
        bader_stats.pop(k, None)

    return bader_stats


def _chargemol_runner(
    path: str, atomic_densities_path: str | None = None
) -> ChargemolSchema:
    """
    Runs a Chargemol (i.e. DDEC6 + CM5) analysis using the VASP output files in
    the given path. This function requires that the chargemol executable, given
    by the name `Chargemol_09_26_2017_linux_parallel`,
    `Chargemol_09_26_2017_linux_serial`, or `chargemol` is in the system PATH
    environment variable. See https://sourceforge.net/projects/ddec/files for
    the Chargemol code.

    Parameters
    ----------
    path
        The path where the VASP output files are located. Must include CHGCAR,
        AECCAR0, AECCAR2, and POTCAR files. These files can be gzip'd or not --
        it doesn't matter. If None, the current working directory is used.
    atomic_densities_path
        The path where the reference atomic densities are located for Chargemol.
        If None, we assume that this directory is defined in an environment
        variable named DDEC6_ATOMIC_DENSITIES_DIR. See the Chargemol
        documentation for more information.

    Returns
    -------
    ChargemolSchema
        Dictionary containing the Chargemol analysis summary
    """

    # Make sure files are present
    relevant_files = ["AECCAR0", "AECCAR2", "CHGCAR", "POTCAR"]
    for f in relevant_files:
        if not Path(path, f).exists() and not Path(path, f"{f}.gz").exists():
            msg = f"Could not find {f} in {path}."
            raise FileNotFoundError(msg)

    # Check environment variable
    if atomic_densities_path is None and "DDEC6_ATOMIC_DENSITIES_DIR" not in os.environ:
        msg = "DDEC6_ATOMIC_DENSITIES_DIR environment variable not defined."
        raise OSError(msg)

    # Run Chargemol analysis
    return ChargemolAnalysis(path=path, atomic_densities_path=atomic_densities_path)
