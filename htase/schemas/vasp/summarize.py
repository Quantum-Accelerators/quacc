from atomate2.vasp.schemas.task import TaskDocument
import os
from ase.io.jsonio import encode
from htase.util.calc import cache_calc
from monty.json import jsanitize


def get_results(atoms=None, dir_path=None, **kwargs):
    """

    Args:
        atoms (ase.Atoms): ASE Atoms object to store in {"atoms": atoms}.
            Defaults to None.
        dir_path (str): Path to VASP outputs
            Defaults to None (current working directory)

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
    )
    for unused_prop in unused_props:
        results.pop(unused_prop, None)

    if atoms:

        # Stores calculator results in the atoms.info flag and moves
        # final magmoms to initial (necessary for sequential jobs)
        # Note: because Atoms objects are mutable, this change will
        # carry through even though we do not return it
        atoms = cache_calc(atoms)

        # Store the encoded Atoms object for later use
        results["atoms"] = encode(atoms)

        results["atoms_info"] = {}
        for key, val in atoms.info.items():
            results["atoms_info"][key] = jsanitize(val)

    return results
