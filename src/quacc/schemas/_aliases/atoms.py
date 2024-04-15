"""Aliases for type hinting `quacc.schemas.atoms`"""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.schemas._aliases.emmet import MoleculeMetadata, StructureMetadata

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from pymatgen.core.structure import Molecule, Structure


class AtomsSchema(StructureMetadata, MoleculeMetadata):
    """Type hint associated with [quacc.schemas.atoms.atoms_to_metadata][]"""

    atoms: Atoms
    structure: Structure  # if atoms.pbc.any()
    molecule: Molecule  # if not atoms.pbc.any()
