"""Schemas for storing metadata about Atoms objects."""

from __future__ import annotations

from typing import TYPE_CHECKING

from emmet.core.structure import MoleculeMetadata, StructureMetadata
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.atoms.core import copy_atoms

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import AtomsSchema


def atoms_to_metadata(
    atoms: Atoms, additional_fields: dict[str, Any] | None = None
) -> AtomsSchema:
    """
    Convert an ASE Atoms object to a dict suitable for storage in MongoDB.

    Parameters
    ----------
    atoms
        ASE Atoms object to store in {"atoms": atoms}
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
    atoms.calc = None

    # Strip the dummy atoms, if present
    if "X" in atoms.get_chemical_symbols():
        del atoms[[atom.index for atom in atoms if atom.symbol == "X"]]

    # Get Pymatgen Structure/Molecule metadata
    if atoms.pbc.any():
        structure = AseAtomsAdaptor().get_structure(atoms)
        structure_metadata = StructureMetadata().from_structure(structure).model_dump()
        results["structure_metadata"] = structure_metadata
    else:
        mol = AseAtomsAdaptor().get_molecule(atoms, charge_spin_check=False)
        molecule_metadata = MoleculeMetadata().from_molecule(mol).model_dump()
        for key in ["charge", "spin_multiplicity", "nelectrons"]:
            del molecule_metadata[key]
        results["molecule_metadata"] = molecule_metadata

    # Store Atoms object
    results["atoms"] = atoms

    return results | additional_fields
