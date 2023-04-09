"""
Schemas for storing metadata about Atoms objects
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Dict

import numpy as np
from ase.atoms import Atom, Atoms
from emmet.core.structure import MoleculeMetadata, StructureMetadata
from monty.json import jsanitize
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.util.atoms import copy_atoms


def atoms_to_metadata(
    atoms: Atoms, get_metadata: bool = True, strip_info: bool = False, store_pmg=True
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
        Note that this data will be stored in {"atoms_info":atoms.info} regardless
    store_pmg
        Whether to store the Pymatgen Structure/Molecule object in {"structure": Structure}
        or {"molecule": Molecule}, respectively.

    Returns
    -------
    Dict
        Dictionary of tabulated atoms object data
    """

    atoms = copy_atoms(atoms)
    results = {}

    # Get Atoms metadata, if requested. Atomate2 already has built-in tools for
    # generating pymatgen Structure/Molecule metadata, so we'll just use that.
    if get_metadata:
        if atoms.pbc.any():
            struct = AseAtomsAdaptor().get_structure(atoms)
            metadata = StructureMetadata().from_structure(struct).dict()
            if store_pmg:
                results["structure"] = struct
        else:
            mol = AseAtomsAdaptor().get_molecule(atoms, charge_spin_check=False)
            metadata = MoleculeMetadata().from_molecule(mol).dict()
            if store_pmg:
                results["molecule"] = mol
        metadata = _quacc_sanitize(metadata)
    else:
        metadata = {}

    # Copy the info flags as a separate entry in the DB for easy querying
    results["atoms_info"] = _quacc_sanitize(atoms.info)

    # Strip info if requested
    if strip_info:
        atoms_no_info = copy_atoms(atoms)
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
    elif isinstance(obj, datetime):
        obj = str(datetime)
    else:
        obj = jsanitize(obj)
    return obj
