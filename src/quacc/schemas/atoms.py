"""Schemas for storing metadata about Atoms objects"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atom, Atoms
from emmet.core.structure import MoleculeMetadata, StructureMetadata
from monty.json import jsanitize
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.utils.dicts import sort_dict

if TYPE_CHECKING:
    from typing import Any, TypedDict

    class AtomsSchema(TypedDict, total=False):
        """Schema for `atoms_to_metadata`"""

        atoms: Atoms
        atoms_info: dict[str, Any]  # from atoms.info
        structure: StructureMetadata  # if atoms.pbc.any()
        molecule: MoleculeMetadata  # if not atoms.pbc.any()


def atoms_to_metadata(
    atoms: Atoms,
    charge_and_multiplicity: tuple[int, int] | None = None,
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

    Returns
    -------
    AtomsSchema
        Dictionary of metadata about the Atoms object.
    """

    additional_fields = additional_fields or {}
    results = {}

    # Set any charge or multiplicity keys. This is a bit of a hack
    # in order to ensure that the Structure or Molecule metadata
    # gets the right charge and multiplicity.
    if charge_and_multiplicity:
        atoms.charge = charge_and_multiplicity[0]
        atoms.spin_multiplicity = charge_and_multiplicity[1]

    # Strip the dummy atoms, if present
    del atoms[[atom.index for atom in atoms if atom.symbol == "X"]]

    # Get Atoms metadata
    results["atoms"] = atoms

    if atoms.pbc.any():
        struct = AseAtomsAdaptor().get_structure(atoms)
        metadata = StructureMetadata().from_structure(struct).dict()
        results["structure"] = struct
    else:
        mol = AseAtomsAdaptor().get_molecule(atoms, charge_spin_check=False)
        metadata = MoleculeMetadata().from_molecule(mol).dict()
        results["molecule"] = mol

    # Copy the info flags as a separate entry in the DB for easy querying
    results["atoms_info"] = _quacc_sanitize(atoms.info)

    # Combine the metadata and results dictionaries
    atoms_doc = metadata | results

    return sort_dict(atoms_doc)


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
