import os
from typing import Any, Dict
from ase.atoms import Atoms
from atomate2.vasp.schemas.task import TaskDocument
from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_
from quacc.util.json import jsanitize


def summarize_run(
    atoms: Atoms,
    dir_path: None | str = None,
    prep_next_run: bool = True,
    check_convergence: bool = True,
    **taskdoc_kwargs
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

    # Remove some key/vals we don't actually ever use
    unused_props = (
        "icsd_id",
        "author",
        "calcs_reversed",
        "transformations",
        "entry",
        "task_label",
        "tags",
        "structure",
    )
    for unused_prop in unused_props:
        results.pop(unused_prop, None)
    if results.get("included_objects", None) is None:
        results.pop("included_objects", None)
    if results.get("vasp_objects", {}) == {}:
        results.pop("vasp_objects", None)

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
    task_doc = jsanitize(results_full)

    return task_doc
