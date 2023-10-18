"""Aliases for type hinting `quacc.schemas.atoms`"""
from typing import Any

from ase import Atoms
from pymatgen.core import Molecule, Structure

from quacc.schemas._aliases.emmet import MoleculeMetadata, StructureMetadata


class AtomsSchema(StructureMetadata, MoleculeMetadata):
    """
    Type hint associated with `quacc.schemas.atoms.atoms_to_metadata`
    """

    atoms: Atoms
    atoms_info: dict[str, Any]  # from atoms.info
    structure: Structure  # if atoms.pbc.any()
    molecule: Molecule  # if not atoms.pbc.any()
