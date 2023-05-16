"""
Schemas for molecular DFT codes parsed by cclib
"""
from __future__ import annotations

import os
from typing import Any, Dict, List

from ase import Atoms
from atomate2.common.schemas.cclib import TaskDocument

from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_


def summarize_run(
    atoms: Atoms,
    logfile_extensions: str | List[str],
    dir_path: str = None,
    check_convergence: bool = True,
    transition_state: bool = False,
    prep_next_run: bool = True,
    additional_fields: Dict[str, Any] = None,
) -> Dict[str, Any]:
    """
    Get tabulated results from a molecular DFT run and store them in a database-friendly format.
    This is meant to be a general parser built on top of cclib.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    logfile_extensions
        Possible extensions of the log file (e.g. ".log", ".out", ".txt", ".chk"). Note that
        only a partial match is needed. For instance, `.log` will match `.log.gz` and `.log.1.gz`.
        If multiple files with this extension are found, the one with the most recent change time
        will be used. For an exact match only, put in the full file name.
    dir_path
        The path to the folder containing the calculation outputs. A value of None specifies the
        current working directory.
    check_convergence
         Whether to throw an error if convergence is not reached.
    transition_state
        Whether the calculation is a transition state (used for convergence check).
    prep_next_run
        Whether the Atoms object storeed in {"atoms": atoms} should be prepared for the next run.
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    additional_fields
        Additional fields to add to the task document.

    Returns
    -------
    Dict
        Dictionary of tabulated inputs/results
    """
    # Make sure there is a calculator with results
    if not atoms.calc:
        raise ValueError("ASE Atoms object has no attached calculator.")
    if not atoms.calc.results:
        raise ValueError("ASE Atoms object's calculator has no results.")

    additional_fields = additional_fields or {}
    dir_path = dir_path or os.getcwd()

    # Fortunately, there is already a cclib parser in Atomate2
    results = TaskDocument.from_logfile(dir_path, logfile_extensions).dict()
    uri = results["dir_name"]
    results["nid"] = uri.split(":")[0]
    results["dir_name"] = ":".join(uri.split(":")[1:])
    results["builder_meta"]["build_date"] = str(results["builder_meta"]["build_date"])

    # Check convergence if requested
    if check_convergence:
        # If it's an opt+freq job, we will just check convergence on the
        # frequency step. This is because sometimes the frequency job will
        # yield all positive modes, but the optimization tolerances may
        # not be entirely met. See https://gaussian.com/faq3.
        vibfreqs = results["attributes"].get("vibfreqs")
        if vibfreqs:
            n_imag = sum(vibfreq < 0 for vibfreq in vibfreqs)
            if n_imag >= 2:
                raise ValueError(f"Too many imaginary modes: {n_imag}")
            if n_imag == 1 and not transition_state:
                raise ValueError("One imaginary mode, but transition_state = False.")
            if n_imag == 0 and transition_state:
                raise ValueError("No imaginary modes, but transition_state = True.")
        elif results["attributes"].get("optdone") is False:
            raise ValueError("Optimization not complete.")

    # Remove some key/vals we don't actually ever use
    unused_props = (
        "task_label",
        "tags",
    )
    for unused_prop in unused_props:
        results.pop(unused_prop, None)

    # Get the calculator inputs
    inputs = {"parameters": atoms.calc.parameters}

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # Get tabulated properties of the structure itself
    atoms_db = atoms_to_metadata(atoms)

    # Create a dictionary of the inputs/outputs
    task_doc = {**atoms_db, **inputs, **results, **additional_fields}

    return task_doc
