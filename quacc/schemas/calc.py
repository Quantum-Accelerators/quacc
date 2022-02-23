from typing import Any, Dict

from ase.atoms import Atoms

from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_


def summarize_run(atoms: Atoms,
                  prep_next_run: bool = True,
                  additional_fields: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Get tabulated results from an Atoms object and calculator and store them in a database-friendly format.
    This is meant to be compatible with all calculator types.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared for the next run
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

    if additional_fields is None:
        additional_fields = {}

    # Fetch all tabulated results from the attached calculator
    results = {"results": atoms.calc.results}

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
