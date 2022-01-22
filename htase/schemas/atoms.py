from atomate2.common.schemas.structure import StructureMetadata
from atomate2.common.schemas.molecule import MoleculeMetadata
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.jsonio import encode
from monty.json import jsanitize
import numpy as np
from copy import deepcopy


def atoms_to_db(atoms, get_metadata=True, strip_info=False):

    """
    Convert an ASE Atoms object to a dict suitable for storage in MongoDB.

    Args:
        atoms (ase.Atoms): ASE Atoms object to store in {"atoms": atoms}.
        get_metadata (bool): Whether to store atoms metadata in the returned dict.
            Defaults to True.
        strip_info (bool): Whether to strip the data from atoms.info in the returned {"atoms":atoms}.
            Note that this data will be stored in {"atoms_info":atoms.info} regardless.
            Defaults to False.

    Returns:
        Dict: dictionary of tabulated atoms object data

    """

    atoms = deepcopy(atoms)
    atoms_no_info = deepcopy(atoms)
    atoms_no_info.info = {}
    results = {}

    # Get Atoms metadata, if requested. Atomate2 already has built-in tools for
    # generating pymatgen Structure/Molecule metadata, so we'll just use that.
    if get_metadata:
        if np.all(atoms.pbc == False):
            mol = AseAtomsAdaptor().get_molecule(atoms)
            metadata = MoleculeMetadata().from_molecule(mol).dict()
        else:
            struct = AseAtomsAdaptor().get_structure(atoms)
            metadata = StructureMetadata().from_structure(struct).dict()
    else:
        metadata = {}

    # Copy the info flags as a separate entry in the DB for easy querying
    results["atoms_info"] = {}
    for key, val in atoms.info.items():
        results["atoms_info"][key] = val

    # Store the encoded Atoms object
    if strip_info:
        results["atoms"] = encode(atoms_no_info)
    else:
        results["atoms"] = encode(atoms)

    # Combine the metadata and results dictionaries
    results_full = {**metadata, **results}

    # Make sure it's all JSON serializable
    results_full = jsanitize(results_full)

    return results_full
