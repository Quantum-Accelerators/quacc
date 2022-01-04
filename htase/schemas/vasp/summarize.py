from atomate2.vasp.schemas.task import TaskDocument
import os
from ase.io.jsonio import encode
from htase.util.calc import cache_calc


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

    if atoms:

        # Store the encoded Atoms object for later use
        results["atoms"] = encode(atoms)

        # Store additional information if present in atoms.info
        for k, v in atoms.info.items():
            results[k] = v

        # Stores calculator results in the atoms.info flag and moves
        # final magmoms to initial (necessary for sequential jobs)
        atoms = cache_calc(atoms)

    return results
