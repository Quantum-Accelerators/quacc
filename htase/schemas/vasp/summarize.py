from atomate2.vasp.schemas.task import TaskDocument
import os
from htase.schemas.common.atoms import atoms_to_db
from htase.util.atoms import prep_next_run as prep_next_run_func


def get_results(atoms, dir_path=None, prep_next_run=True, **taskdoc_kwargs):
    """

    Args:
        atoms (ase.Atoms): ASE Atoms object following a calculation.
        dir_path (str): Path to VASP outputs
            Defaults to None (current working directory)
        prep_next_run (bool): Whether the Atoms object storeed in {"atoms": atoms} should be prepared
            for the next run. This clears out any attached calculator and moves the final magmoms to the
            initial magmoms.
            Defauls to True.
        **taskdoc_kwargs: Additional keyword arguments to pass to TaskDocument.from_directory()

    Returns:
        results (dict): dictionary of tabulated results

    """

    if dir_path is None:
        dir_path = os.getcwd()

    # Fetch all tabulated results from VASP outputs files
    # Fortunately, Atomate2 already has a handy function for this
    results = TaskDocument.from_directory(dir_path, **taskdoc_kwargs).dict()

    # Remove some key/vals we don't actually ever use
    unused_props = (
        "icsd_id",
        "author",
        "calcs_reversed",
        "transformations",
        "state",
        "entry",
        "task_label",
        "tags",
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
        atoms = prep_next_run_func(atoms)

    # We use get_metadata=False because the TaskDocument already
    # makes the structure metadata for us
    atoms_db = atoms_to_db(atoms, get_metadata=False)

    results_full = {**results, **atoms_db}

    return results_full
