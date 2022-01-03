from atomate2.vasp.schemas.task import TaskDocument
import os


def get_results(dir_path=None, atoms=None, **kwargs):
    """

    Args:
        dir_path (str): Path to VASP outputs
            Defaults to None (current working directory)
        atoms (ase.Atoms): ASE Atoms object to store in {"atoms": atoms}
            Defaults to None (not stored)
        **kwargs: additional keyword arguments to pass to TaskDocument

    Returns:
        results (dict): dictionary of tabulated results

    """

    if dir_path is None:
        dir_path = os.getcwd()
    results = TaskDocument.from_directory(dir_path, **kwargs)
    if atoms is not None:
        results["atoms"] = atoms

    return results
