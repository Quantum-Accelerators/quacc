"""Aliases for type hinting `quacc.schemas.atoms`"""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, TypedDict

    from ase import Atoms
    from pymatgen.core import Molecule, Structure

    from quacc.schemas._aliases.emmet import MoleculeMetadata, StructureMetadata

    class AtomsSchema(TypedDict, StructureMetadata, MoleculeMetadata, total=False):
        """
        Type hint associated with `quacc.schemas.atoms.atoms_to_metadata`
        """

        atoms: Atoms
        atoms_info: dict[str, Any]  # from atoms.info
        structure: Structure  # if atoms.pbc.any()
        molecule: Molecule  # if not atoms.pbc.any()
