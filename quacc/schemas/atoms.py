from copy import deepcopy
from typing import Any, Dict

import numpy as np
from ase.atoms import Atom, Atoms
from atomate2.common.schemas.molecule import MoleculeMetadata
from atomate2.common.schemas.structure import StructureMetadata
from monty.json import jsanitize
from pymatgen.io.ase import AseAtomsAdaptor


def atoms_to_metadata(
    atoms: Atoms, get_metadata: bool = True, strip_info: bool = False
) -> Dict[str, Any]:

    """
    Convert an ASE Atoms object to a dict suitable for storage in MongoDB.

    Parameters
    ----------
    atoms
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
    results["atoms_info"] = _quacc_sanitize(atoms.info)

    # Strip info if requested
    if strip_info:
        atoms_no_info = deepcopy(atoms)
        atoms_no_info.info = {}
        results["atoms"] = atoms_no_info
    else:
        results["atoms"] = atoms

    # Combine the metadata and results dictionaries
    atoms_doc = {**metadata, **results}

    return atoms_doc


def _quacc_sanitize(obj: Any) -> Any:
    """
    Sanitizes an object for storage in MongoDB.

    This is an analogue of monty's jsanitize function but meant to serialize
    Atom/Atoms objects as well.

    Parameters
    ----------
    obj
        Object to sanitize

    Returns
    -------
    Any
        Sanitized object
    """
    if isinstance(obj, (Atom, Atoms)):
        obj = atoms_to_metadata(obj)
    elif isinstance(obj, (list, tuple, np.ndarray)):
        obj = [_quacc_sanitize(i) for i in obj]
    elif isinstance(obj, dict):
        obj = {k.__str__(): _quacc_sanitize(v) for k, v in obj.items()}
    else:
        obj = jsanitize(obj)
    return obj
