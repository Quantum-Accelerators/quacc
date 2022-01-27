import os
from typing import Dict, Optional
from ase.atoms import Atoms
from cclib.io import ccopen
from monty.json import jsanitize
from quacc.schemas.atoms import atoms_to_db
from quacc.util.atoms import prep_next_run as prep_next_run_


def summarize_run(
    atoms: Atoms, output_file: Optional[str] = None, prep_next_run: bool = True
) -> Dict:
    """
    Get tabulated results from a molecular DFT run and store them in a database-friendly format.
    This is meant to be a general parser built on top of cclib.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    output_file
        Path to the main output file. If None, common output files will be searched for automatically
        in the working directory, and the one with the most recent modification time will be used.
    prep_next_run
        Whether the Atoms object storeed in {"atoms": atoms} should be prepared for the next run.
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.

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

    # Try to find the most recent output file
    if output_file is None:
        logfile_list = [".log", ".out", ".chk"]  # needs to be updated
        mod_time = 0
        for f in os.listdir(os.getcwd()):
            if f.split(".gz")[0] in logfile_list:
                if output_file:
                    if os.path.getmtime(f) > mod_time:
                        output_file = f
                        mod_time = os.path.getmtime(output_file)
                else:
                    output_file = os.path.getmtime(f)
                    mod_time = os.path.getmtime(output_file)

    # Fetch all tabulated results from the attached calculator
    if not output_file or not os.path.exists(output_file):
        raise FileNotFoundError(f"Could not find {output_file}")
    output_io = ccopen(output_file)
    if output_io is None:
        raise ValueError(f"cclib could not parse {output_file}")

    output = output_io.parse()
    results = {"attributes": output.getattributes()}
    results["attributes"].pop("metadata")
    metadata = {"metadata": output.metadata}

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
    results_full = {**atoms_db, **inputs, **metadata, **results}

    # Make sure it's all JSON serializable
    task_doc = jsanitize(results_full)

    return task_doc
