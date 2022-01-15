from atomate2.vasp.schemas.task import TaskDocument
import os
from ase.io.jsonio import encode
from htase.util.calc import cache_calc
from monty.json import jsanitize
from copy import deepcopy


def get_results(atoms=None, dir_path=None, tags=None, **kwargs):
    """

    Args:
        atoms (ase.Atoms): ASE Atoms object to store in {"atoms": atoms}.
            Defaults to None.
        dir_path (str): Path to VASP outputs
            Defaults to None (current working directory)
        tags (List[str]): List of tags to store in {"tags": tags}. Useful
            for storing job metadata.

    Returns:
        results (dict): dictionary of tabulated results

    """

    if dir_path is None:
        dir_path = os.getcwd()

    # Fetch all tabulated results from VASP outputs files
    # Fortunately, Atomate2 already has a handy function for this
    results = TaskDocument.from_directory(dir_path, **kwargs).dict()

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

    # Clean up
    if results["included_objects"] is None:
        del results["included_objects"]
    if results["vasp_objects"] == {}:
        del results["vasp_objects"]

    if atoms:

        # Stores calculator results in the atoms.info flag and moves
        # final magmoms to initial (necessary for sequential jobs)
        # Note: because Atoms objects are mutable, this change will
        # carry through even though we do not return the Atoms object
        atoms = cache_calc(atoms)

        # Store the info flags separately
        atoms_noinfo = deepcopy(atoms)
        atoms_noinfo.info = {}

        # Store the info flags separately
        results["info"] = {}
        for key, val in atoms.info.items():
            results["info"][key] = jsanitize(val)

        # Store any tags
        if tags:
            results["tags"] = tags

        # Store the encoded Atoms object (without .info)
        results["atoms"] = encode(atoms_noinfo)

    return results
