import os
from typing import Any, Dict
from ase.atoms import Atoms
from atomate2.vasp.schemas.task import TaskDocument
from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_
from quacc.util.json import clean


def summarize_run(
    atoms: Atoms,
    dir_path: None | str = None,
    prep_next_run: bool = True,
    check_convergence: bool = True,
    compact: bool = True,
    remove_empties: bool = True,
    **taskdoc_kwargs,
) -> Dict[str, Any]:
    """
    Get tabulated results from a VASP run and store them in a database-friendly format.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    dir_path
        Path to VASP outputs. A value of None specifies the current working directory
    prep_next_run
        Whether the Atoms object storeed in {"atoms": atoms} should be prepared for the next run.
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    check_convergence
        Whether to throw an error if convergence is not reached.
    compact
        Whether to use a compact representation of the TaskDocument. This does not store the
        inputs/outputs from intermediate Custodian runs. It also removes duplicate key-value
        pairs.
    remove_empties
        Whether to remove None values and empty lists/dicts from the TaskDocument.
    **taskdoc_kwargs: Additional keyword arguments to pass to TaskDocument.from_directory()

    Returns
    -------
    Dict
        Dictionary of tabulated inputs/results
    """

    if dir_path is None:
        dir_path = os.getcwd()

    # Fetch all tabulated results from VASP outputs files
    # Fortunately, Atomate2 already has a handy function for this
    results = TaskDocument.from_directory(dir_path, **taskdoc_kwargs).dict()

    # Check for calculation convergence
    if check_convergence and results["state"] != "successful":
        raise RuntimeError(
            "VASP calculation did not converge. Will not store task data."
        )

    if compact:
        # Replace the InputSummary and OutputSummary with the full
        # input and output details from calcs_reversed
        if results.get("calcs_reversed", None):
            final_run = results["calcs_reversed"][-1]
            results["input"] = final_run["input"]
            results["output"] = final_run["output"]

            # Store a few additional properties
            results["vasp_version"] = final_run["vasp_version"]
            results["task_type"] = final_run["task_type"]
            results["run_type"] = final_run["run_type"]

            # Then dump the calcs_reversed
            results.pop("calcs_reversed")

        # Remove structure because it's already in the outputs
        results.pop("structure", None)

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # We use get_metadata=False because the TaskDocument already
    # makes the structure metadata for us
    atoms_db = atoms_to_metadata(atoms, get_metadata=False)

    results_full = {**results, **atoms_db}

    # Make sure it's all JSON serializable
    task_doc = clean(results_full)

    if remove_empties:
        task_doc = _remove_empties(task_doc)

    return task_doc


def _remove_empties(d: Dict[str, Any]) -> Dict[str, Any]:
    """
    For a given dictionary, recursively remove all items that are None
    or are empty lists/dicts.

    Parameters
    ----------
    d
        Dictionary to clean

    Returns
    -------
    Dict
        Cleaned dictionary
    """

    if isinstance(d, dict):
        return {
            k: _remove_empties(v)
            for k, v in d.items()
            if v != None and v != [] and v != {}
        }
    if isinstance(d, list):
        return [_remove_empties(v) for v in d]
    return d
