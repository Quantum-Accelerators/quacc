from ase.atoms import Atoms
from ase.io.jsonio import encode
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from copy import deepcopy
import hashlib
import os

# NOTES:
# - Anytime an Atoms object is converted to a pmg structure, make sure
# to reattach any .info flags to the Atoms object, e.g. via `new_atoms.info = atoms.info.copy()``.
# Note that atoms.info is mutable, so copy it!
# - All major functions should take in Atoms by default and reutrn Atoms
# by default. Pymatgen structures can be returned with an optional kwarg.
# - If you modify the properties of an input Atoms object in any way, make sure to do so
# on a deepcopy because Atoms objects are mutable.
# - If you are going to store an Atoms/Atom object in the atoms.info dictionary, do so using
# atoms_to_db(atoms) so that it can be properly serialized.


def prep_next_run(
    atoms,
    assign_id=True,
    move_magmoms=True,
    store_results=False,
):
    """
    Prepares the Atoms object for a new run.

    Depending on the arguments, this function will:
        - Move the converged magnetic moments to the initial magnetic moments.
        - Assign a unique ID to the Atoms object in atoms.info["_id"]. Any existing IDs will
        be moved to atoms.info["_old_ids"].
        - Store the calculator results in atoms.info["results"] for later retrieval.
        This makes it so the calculator results are not lost between
        serialize/deserialize cycles, if desired. Each one will be stored in
        atoms.info["results"] = {"calc0": {}, "calc1": {}, ...} with higher numbers
        being the most recent.

    In all cases, the calculator will be reset so new jobs can be run.

    Args:
        atoms (ase.Atoms): Atoms object
        assign_id (bool): Whether to assign a unique ID to the Atoms object in atoms.info["_id"].
            Any existing IDs will be moved to atoms.info["_old_ids"].
            move_magmoms (bool): If True, move atoms.get_magnetic_moments() to
            atoms.get_initial_magnetic_moments()
            Defult: True.
        store_results (bool): If True, store calculator results in atoms.info["results"].
            This makes it so the calculator results are not lost between serialize/deserialize cycles,
            if desired. Each one will be stored in atoms.info["results"] = {"calc0": {}, "calc1": {}, ...}
            with higher numbers being the most recent.
            Default: False.

    Returns:
        atoms (ase.Atoms): Atoms object with calculator results attached in atoms.info["results"]
    """
    atoms = deepcopy(atoms)

    if hasattr(atoms, "calc") and getattr(atoms.calc, "results", None) is not None:

        if store_results:
            # Dump calculator results into the .info tag
            atoms.calc.results["rundir"] = os.getcwd()
            if atoms.info.get("results", None) is None:
                prior_calcs = 0
                atoms.info["results"] = {}
            else:
                prior_calcs = len(atoms.info["results"])

            atoms.info["results"][f"calc{prior_calcs}"] = atoms.calc.results

        # Move converged magmoms to initial magmoms
        # If none were present, then initial magmoms should be set to 0's
        # because a spin-unpolarized calculation was carried out
        if move_magmoms:
            atoms.set_initial_magnetic_moments(
                atoms.calc.results.get("magmoms", [0.0] * len(atoms))
            )

    # Clear off the calculator so we can run a new job
    atoms.calc = None

    # Give the Atoms object a unique ID. This will be helpful for querying later.
    # Also store any old IDs somewhere else for future reference.
    if assign_id:
        if atoms.info.get("_id", None) is not None:
            if atoms.info.get("_old_ids") is None:
                atoms.info["_old_ids"] = []
            atoms.info["_old_ids"].append(atoms.info["_id"])
        atoms.info["_id"] = get_atoms_id(atoms)

    return atoms


def get_atoms_id(atoms):
    """
    Returns a unique ID for the Atoms object. Note: The .info dict
    and calculator is excluded from the hash generation.

    Args:
        atoms (ase.Atoms): Atoms object

    Returns:
        md5hash (str): MD5 hash of the Atoms object
    """

    atoms = deepcopy(atoms)
    atoms.info = {}
    atoms.calc = None
    encoded_atoms = encode(atoms)
    # This is a hack to avoid int32/int64 and float32/float64 differences
    # between machines.
    encoded_atoms = (
        encoded_atoms.replace("int64", "int")
        .replace("int32", "int")
        .replace("float64", "float")
        .replace("float32", "float")
    )
    md5hash = hashlib.md5(encoded_atoms.encode("utf-8")).hexdigest()
    return md5hash


def check_is_metal(atoms):
    """
    Checks if a structure is a likely metal.

    Args:
        atoms (ase.Atoms): ASE atoms object

    Returns:
        is_metal (bool): True if the structure is likely a metal; False otherwise
    """
    if type(atoms) is Atoms:
        if np.all(atoms.pbc) == False:
            struct = AseAtomsAdaptor.get_molecule(atoms)
        else:
            struct = AseAtomsAdaptor.get_structure(atoms)
    else:
        struct = atoms
    is_metal = all(k.is_metal for k in struct.composition.keys())
    return is_metal


def get_highest_block(atoms):
    """
    Get the highest block (e.g. p-block, d-block f-block) of a structure

    Args:
        atoms (ase.Atoms): ASE atoms object

    Returns:
        highest_block (str): highest block of the structure
    """
    if type(atoms) is Atoms:
        if np.all(atoms.pbc) == False:
            struct = AseAtomsAdaptor.get_molecule(atoms)
        else:
            struct = AseAtomsAdaptor.get_structure(atoms)
    else:
        struct = atoms
    blocks = [site.specie.block for site in struct]
    if "f" in blocks:
        max_block = "f"
    elif "d" in blocks:
        max_block = "d"
    elif "p" in blocks:
        max_block = "p"
    else:
        max_block = "s"

    return max_block
