"""Schemas for storing metadata about Atoms objects."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atom, Atoms
from emmet.core.structure import MoleculeMetadata, StructureMetadata
from monty.json import jsanitize
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.atoms.core import copy_atoms
from quacc.utils.dicts import remove_dict_nones, sort_dict

if TYPE_CHECKING:
    from typing import Any

    from quacc.schemas._aliases.atoms import AtomsSchema


def atoms_to_metadata(
    atoms: Atoms,
    charge_and_multiplicity: tuple[int, int] | None = None,
    get_metadata: bool = True,
    store_pmg: bool = True,
    additional_fields: dict[str, Any] | None = None,
) -> AtomsSchema:
    """
    Convert an ASE Atoms object to a dict suitable for storage in MongoDB.

    Parameters
    ----------
    atoms
        ASE Atoms object to store in {"atoms": atoms}
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule
        metadata.
    get_metadata
        Whether to store atoms metadata in the returned dict.
    store_pmg
        Whether to store the Pymatgen Structure/Molecule object in {"structure":
        Structure} or {"molecule": Molecule}, respectively.
    additional_fields
        Additional fields to add to the document.

    Returns
    -------
    AtomsSchema
        Dict of metadata about the Atoms object.
    """

    additional_fields = additional_fields or {}
    atoms = copy_atoms(atoms)
    results = {}

    # Get any charge or multiplicity keys
    if charge_and_multiplicity:
        atoms.charge = charge_and_multiplicity[0]
        atoms.spin_multiplicity = charge_and_multiplicity[1]

    # Strip the dummy atoms, if present
    del atoms[[atom.index for atom in atoms if atom.symbol == "X"]]

    # Get Atoms metadata, if requested. emmet already has built-in tools for
    # generating pymatgen Structure/Molecule metadata, so we'll just use that.
    if get_metadata:
        if atoms.pbc.any():
            struct = AseAtomsAdaptor().get_structure(atoms)
            metadata = StructureMetadata().from_structure(struct).model_dump()
            if store_pmg:
                results["structure"] = struct
        else:
            mol = AseAtomsAdaptor().get_molecule(atoms, charge_spin_check=False)
            metadata = MoleculeMetadata().from_molecule(mol).model_dump()
            if store_pmg:
                results["molecule"] = mol
    else:
        metadata = {}

    # Copy the info flags as a separate entry in the DB for easy querying
    results["atoms_info"] = _quacc_sanitize(atoms.info)

    # Store Atoms object
    results["atoms"] = atoms

    # Combine the metadata and results dictionaries
    atoms_doc_unsorted = metadata | results | additional_fields

    return sort_dict(remove_dict_nones(atoms_doc_unsorted))


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
