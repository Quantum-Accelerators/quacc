from monty.json import jsanitize
from htase.schemas.atoms import atoms_to_db
from htase.util.atoms import prep_next_run as prep_next_run_


def summarize_run(atoms, prep_next_run=True):
    """
    Get tabulated results from an Atoms object and calculator and store them in a database-friendly format.
    This is meant to be compatible with all calculator types.

    Args:
        atoms (ase.Atoms): ASE Atoms object following a calculation.
        prep_next_run (bool): Whether the Atoms object storeed in {"atoms": atoms} should be prepared
            for the next run. This clears out any attached calculator and moves the final magmoms to the
            initial magmoms.
            Defauls to True.

    Returns:
        task_doc (dict): dictionary of tabulated inputs/results

    """
    # Make sure there is a calculator with results
    if not atoms.calc:
        raise ValueError("ASE Atoms object has no attached calculator.")
    if not atoms.calc.results:
        raise ValueError("ASE Atoms object's calculator has no results.")

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
    atoms_db = atoms_to_db(atoms)

    # Create a dictionary of the inputs/outputs
    results_full = {**atoms_db, **inputs, **results}

    # Make sure it's all JSON serializable
    task_doc = jsanitize(results_full)

    return task_doc
