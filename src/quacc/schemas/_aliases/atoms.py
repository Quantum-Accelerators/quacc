"""Aliases for type hinting `quacc.schemas.atoms`"""

from __future__ import annotations

from typing import Any

from ase.atoms import Atoms
from pymatgen.core.structure import Molecule, Structure

from quacc.schemas._aliases.emmet import MoleculeMetadata, StructureMetadata


class AtomsSchema(StructureMetadata, MoleculeMetadata):
    """Type hint associated with [quacc.schemas.atoms.atoms_to_metadata][]"""

    atoms: Atoms
    atoms_info: dict[str, Any]  # from atoms.info
    structure: Structure  # if atoms.pbc.any()
    molecule: Molecule  # if not atoms.pbc.any()
