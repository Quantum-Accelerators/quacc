"""Schemas for VASP"""
from __future__ import annotations

import os
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

from emmet.core.tasks import TaskDoc
from maggma.core import Store
from pymatgen.command_line.bader_caller import bader_analysis_from_path
from pymatgen.command_line.chargemol_caller import ChargemolAnalysis

from quacc import SETTINGS
from quacc.runners.prep import prep_next_run as prep_next_run_
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import sort_dict
from quacc.utils.files import copy_decompress
from quacc.wflow.db import results_to_db

if TYPE_CHECKING:
    from typing import TypedDict, TypeVar

    from ase import Atoms

    VaspSchema = TypeVar("VaspSchema", dict)

    class BaderSchema(TypedDict):
        atomic_volume: float
        bader_charge: float
        bader_spin: float
        bader_version: float
        min_dist: list[float]
        partial_charges: list[float]
        spin_moments: list[float]

    class DDECSchema(TypedDict):
        partial_charges: list[float]
        spin_moments: list[float]
        dipoles: list[float]
        bond_order_sums: list[float]
        bond_order_dict: dict
        rsquared_moments: list[float]
        rcubed_moments: list[float]
        rfourth_moments: list[float]


def vasp_summarize_run(
    atoms: Atoms,
    dir_path: str | None = None,
    prep_next_run: bool = True,
    run_bader: bool | None = None,
    check_convergence: bool = True,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> VaspSchema:
    """
    Get tabulated results from a VASP run and store them in a database-friendly
    format.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    dir_path
        Path to VASP outputs. A value of None specifies the current working
        directory
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared
        for the next run. This clears out any attached calculator and moves the
        final magmoms to the initial magmoms.
    run_bader
        Whether a Bader analysis should be performed. Will not run if bader
        executable is not in PATH even if bader is set to True. Defaults to
        VASP_BADER in settings.
    check_convergence
        Whether to throw an error if convergence is not reached.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    VaspSchema
        Dictionary representation of the task document
    """

    additional_fields = additional_fields or {}
    run_bader = SETTINGS.VASP_BADER if run_bader is None else run_bader
    dir_path = dir_path or Path.cwd()
    store = SETTINGS.PRIMARY_STORE if store is None else store

    # Fetch all tabulated results from VASP outputs files Fortunately, emmet
    # already has a handy function for this
    taskdoc = TaskDoc.from_directory(dir_path).dict()

    uri = taskdoc["dir_name"]
    taskdoc["nid"] = uri.split(":")[0]
    taskdoc["dir_name"] = ":".join(uri.split(":")[1:])
    taskdoc["builder_meta"]["build_date"] = str(taskdoc["builder_meta"]["build_date"])

    # Check for calculation convergence
    if check_convergence and taskdoc["state"] != "successful":
        raise ValueError("VASP calculation did not converge. Will not store task data.")

    # Remove unnecessary fields
    for k in [
        "additional_json",
        "author",
        "calcs_reversed",
        "icsd_id",
        "last_updated",
        "structure",  # already in output
        "tags",
        "task_id",
        "task_label",
        "transformations",
        "vasp_objects",
    ]:
        taskdoc.pop(k, None)

    if "output" in taskdoc:
        taskdoc["output"].pop("elph_displaced_structures", None)
        taskdoc["output"].pop("frequency_dependent_dielectric", None)

    # Get Bader analysis
    if run_bader:
        try:
            bader_stats = bader_runner(dir_path)
        except Exception:
            bader_stats = None
            warnings.warn("Bader analysis could not be performed.", UserWarning)

        if bader_stats:
            taskdoc["bader"] = bader_stats

            # Attach bader charges/spins to structure object
            struct = taskdoc["output"]["structure"]
            struct.add_site_property("bader_charge", bader_stats["partial_charges"])
            if "spin_moments" in bader_stats:
                struct.add_site_property("bader_spin", bader_stats["spin_moments"])
            taskdoc["output"]["structure"] = struct

    # Prepares the Atoms object for the next run by moving the final magmoms to
    # initial, clearing the calculator state, and assigning the resulting Atoms
    # object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # We use get_metadata=False and store_pmg=False because the TaskDocument
    # already makes the structure metadata for us
    atoms_db = atoms_to_metadata(atoms, get_metadata=False, store_pmg=False)

    # Make task document
    summary = sort_dict(taskdoc | atoms_db | additional_fields)

    # Store the results
    if store:
        results_to_db(store, taskdoc)

    return summary


def bader_runner(path: str | None = None, scratch_dir: str | None = None) -> dict:
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
    scratch_dir
        The path where the Bader analysis will be run. Defaults to
        SETTINGS.SCRATCH_DIR.

    Returns
    -------
    BaderSchema
        Dictionary containing the Bader analysis summary
    """
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    path = path or Path.cwd()

    # Make sure files are present
    relevant_files = ["AECCAR0", "AECCAR2", "CHGCAR", "POTCAR"]
    for f in relevant_files:
        if not Path(path, f).exists() and not Path(path, f"{f}.gz").exists():
            msg = f"Could not find {f} in {path}."
            raise FileNotFoundError(msg)

    # Run Bader analysis
    with TemporaryDirectory(dir=scratch_dir) as tmpdir:
        copy_decompress(relevant_files, tmpdir)
        bader_stats = bader_analysis_from_path(path)

    # Store the partial charge, which is much more useful than the raw charge
    # and is more intuitive than the charge transferred. An atom with a positive
    # partial charge is cationic, whereas an atom with a negative partial charge
    # is anionic.
    bader_stats["partial_charges"] = [-c for c in bader_stats["charge_transfer"]]

    # Some cleanup of the returned dictionary
    if "magmom" in bader_stats:
        bader_stats["spin_moments"] = bader_stats["magmom"]
    bader_stats.pop("charge", None)
    bader_stats.pop("charge_transfer", None)
    bader_stats.pop("reference_used", None)
    bader_stats.pop("magmom", None)

    return bader_stats


def chargemol_runner(
    path: str | None = None,
    atomic_densities_path: str | None = None,
    scratch_dir: str | None = None,
) -> DDECSchema:
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
    scratch_dir
        The path where the Chargemol analysis will be run.

    Returns
    -------
    DDECSchema
        Dictionary containing the Chargemol analysis summary
    """
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    path = path or Path.cwd()

    # Make sure files are present
    relevant_files = ["AECCAR0", "AECCAR2", "CHGCAR", "POTCAR"]
    for f in relevant_files:
        if not Path(path, f).exists() and not Path(path, f"{f}.gz").exists():
            msg = f"Could not find {f} in {path}."
            raise FileNotFoundError(msg)

    # Check environment variable
    if atomic_densities_path is None and "DDEC6_ATOMIC_DENSITIES_DIR" not in os.environ:
        msg = "DDEC6_ATOMIC_DENSITIES_DIR environment variable not defined."
        raise EnvironmentError(msg)

    # Run Chargemol analysis
    with TemporaryDirectory(dir=scratch_dir) as tmpdir:
        copy_decompress(relevant_files, tmpdir)
        chargemol_stats = ChargemolAnalysis(
            path=path,
            atomic_densities_path=atomic_densities_path,
        )

    return chargemol_stats
