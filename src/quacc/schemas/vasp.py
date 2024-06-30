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

from quacc import get_settings
from quacc.atoms.core import get_final_atoms_from_dynamics
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import finalize_dict, recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.io import Trajectory
    from ase.optimize.optimize import Optimizer
    from maggma.core import Store

    from quacc.schemas._aliases.vasp import (
        BaderSchema,
        ChargemolSchema,
        VaspASEOptSchema,
        VaspSchema,
    )

logger = logging.getLogger(__name__)

_DEFAULT_SETTING = ()


def vasp_summarize_run(
    final_atoms: Atoms,
    directory: str | Path | None = None,
    move_magmoms: bool = True,
    run_bader: bool = _DEFAULT_SETTING,
    run_chargemol: bool = _DEFAULT_SETTING,
    check_convergence: bool = _DEFAULT_SETTING,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    store: Store | None = _DEFAULT_SETTING,
) -> VaspSchema:
    """
    Get tabulated results from a VASP run and store them in a database-friendly format.

    Parameters
    ----------
    final_atoms
        ASE Atoms object following a calculation.
    directory
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
        Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`,

    Returns
    -------
    VaspSchema
        Dictionary representation of the task document
    """
    settings = get_settings()
    run_bader = settings.VASP_BADER if run_bader == _DEFAULT_SETTING else run_bader
    run_chargemol = (
        settings.VASP_CHARGEMOL if run_chargemol == _DEFAULT_SETTING else run_chargemol
    )
    check_convergence = (
        settings.CHECK_CONVERGENCE
        if check_convergence == _DEFAULT_SETTING
        else check_convergence
    )
    directory = Path(directory or final_atoms.calc.directory)
    store = settings.STORE if store == _DEFAULT_SETTING else store
    additional_fields = additional_fields or {}

    # Fetch all tabulated results from VASP outputs files. Fortunately, emmet
    # already has a handy function for this
    vasp_task_model = TaskDoc.from_directory(directory)

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
            f"VASP calculation did not converge. Will not store task data. Refer to {directory}"
        )
    poscar_path = directory / "POSCAR"
    initial_atoms = read(zpath(str(poscar_path)))
    base_task_doc = summarize_run(
        final_atoms, initial_atoms, move_magmoms=move_magmoms, store=None
    )

    if nsteps := len([f for f in os.listdir(directory) if f.startswith("step")]):
        intermediate_vasp_task_docs = {
            "steps": {
                n: TaskDoc.from_directory(directory / f"step{n}").model_dump()
                for n in range(nsteps)
                if (directory / f"step{n}").is_dir()
            }
        }
    else:
        intermediate_vasp_task_docs = {}

    # Get Bader analysis
    if run_bader:
        try:
            bader_results = _bader_runner(directory)
        except Exception:
            bader_results = None
            logging.warning("Bader analysis could not be performed.", exc_info=True)

        if bader_results:
            vasp_task_doc["bader"] = bader_results

    # Get the Chargemol analysis
    if run_chargemol:
        try:
            chargemol_results = _chargemol_runner(directory)
        except Exception:
            chargemol_results = None
            logging.warning("Chargemol analysis could not be performed.", exc_info=True)

        if chargemol_results:
            vasp_task_doc["chargemol"] = chargemol_results

    # Make task document
    unsorted_task_doc = (
        intermediate_vasp_task_docs | vasp_task_doc | base_task_doc | additional_fields
    )
    return finalize_dict(
        unsorted_task_doc, directory, gzip_file=settings.GZIP_FILES, store=store
    )


def summarize_vasp_opt_run(
    optimizer: Optimizer,
    trajectory: Trajectory | list[Atoms] | None = None,
    directory: str | Path | None = None,
    move_magmoms: bool = True,
    run_bader: bool = _DEFAULT_SETTING,
    run_chargemol: bool = _DEFAULT_SETTING,
    check_convergence: bool = _DEFAULT_SETTING,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    store: Store | None = _DEFAULT_SETTING,
) -> VaspASEOptSchema:
    """
    Merges the `vasp_summarize_run` with an `summarize_opt_run`, meant to
    be used for an ASE-based VASP relaxation.

    Parameters
    ----------
    optimizer
        The ASE optimizer object
    trajectory
        ASE Trajectory object or list[Atoms] from reading a trajectory file. If
        None, the trajectory must be found in dyn.trajectory.filename.
    directory
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
        Maggma Store object to store the results in. Defaults to `QuaccSettings.STORE`,
    """
    settings = get_settings()
    store = settings.STORE if store == _DEFAULT_SETTING else store

    final_atoms = get_final_atoms_from_dynamics(optimizer)
    directory = Path(directory or final_atoms.calc.directory)
    opt_run_summary = summarize_opt_run(
        optimizer,
        trajectory=trajectory,
        check_convergence=check_convergence,
        move_magmoms=move_magmoms,
        additional_fields=additional_fields,
        store=None,
    )
    vasp_summary = vasp_summarize_run(
        final_atoms,
        directory=directory,
        move_magmoms=move_magmoms,
        run_bader=run_bader,
        run_chargemol=run_chargemol,
        check_convergence=check_convergence,
        report_mp_corrections=report_mp_corrections,
        additional_fields=additional_fields,
        store=None,
    )
    unsorted_task_doc = recursive_dict_merge(vasp_summary, opt_run_summary)
    return finalize_dict(
        unsorted_task_doc, directory, gzip_file=settings.GZIP_FILES, store=store
    )


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
