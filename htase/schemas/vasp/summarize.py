from atomate2.vasp.schemas.task import TaskDocument
import os
from ase.io.jsonio import encode
from ase.atoms import Atoms


def get_results(dir_path=None, atoms=None, **kwargs):
    """

    Args:
        dir_path (str): Path to VASP outputs
            Defaults to None (current working directory)
        atoms (ase.Atoms): ASE Atoms object to store in {"atoms": atoms}.
        Can be either encoded or decoded.
            Defaults to None (not stored)
        **kwargs: additional keyword arguments to pass to TaskDocument

    Returns:
        results (dict): dictionary of tabulated results

    """

    if dir_path is None:
        dir_path = os.getcwd()
    results = TaskDocument.from_directory(dir_path, **kwargs).dict
    if atoms is not None:
        if isinstance(atoms, Atoms):
            atoms = encode(atoms)
        results["atoms"] = atoms

    return results
