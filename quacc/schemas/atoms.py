from copy import deepcopy
from typing import Any, Dict
import numpy as np
from ase.atoms import Atom, Atoms
from atomate2.common.schemas.structure import StructureMetadata
from atomate2.common.schemas.molecule import MoleculeMetadata
from pymatgen.io.ase import AseAtomsAdaptor
from quacc.util.json import jsanitize


def atoms_to_metadata(
    atoms: Atoms, get_metadata: bool = True, strip_info: bool = False
) -> Dict[str, Any]:

    """
    Convert an ASE Atoms object to a dict suitable for storage in MongoDB.

    Parameters
    ----------
    aomss
        ASE Atoms object to store in {"atoms": atoms}
    get_metadata
        Whether to store atoms metadata in the returned dict.
    strip_info
        Whether to strip the data from atoms.info in the returned {"atoms":.Atoms}.
        Note that this data will be stored in {"atoms_info":atoms.info} regardless.

    Returns
    -------
    Dict
        Dictionary of tabulated atoms object data
    """

    atoms = deepcopy(atoms)
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

        # This is to make sure any Atom/Atoms object in atoms.info is also
        # turned into metadata
        if isinstance(val, (Atom, Atoms)):
            val = atoms_to_metadata(val)
        elif isinstance(val, (list, tuple, np.ndarray)):
            val = [
                atoms_to_metadata(v) if isinstance(v, (Atom, Atoms)) else v for v in val
            ]
        elif isinstance(val, dict):
            val = {
                k.__str__(): atoms_to_metadata(v) if isinstance(v, (Atom, Atoms)) else v
                for k, v in val.items()
            }

        results["atoms_info"][key] = val

    # Strip info if requested
    if strip_info:
        atoms_no_info = deepcopy(atoms)
        atoms_no_info.info = {}
        results["atoms"] = atoms_no_info
    else:
        results["atoms"] = atoms

    # Combine the metadata and results dictionaries
    results_full = {**metadata, **results}

    # Make sure it's all JSON serializable
    results_full = jsanitize(results_full)

    return results_full
