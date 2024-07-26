"""Schemas for storing metadata about Atoms objects."""

from __future__ import annotations

from typing import TYPE_CHECKING

from emmet.core.structure import MoleculeMetadata, StructureMetadata
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.atoms.core import (
    copy_atoms,
    get_charge_attribute,
    get_spin_multiplicity_attribute,
)

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import AtomsSchema


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
    atoms.calc = None

    # Set any charge or multiplicity keys
    if not atoms.pbc.any():
        if charge_and_multiplicity:
            charge = charge_and_multiplicity[0]
            spin_multiplicity = charge_and_multiplicity[1]
        else:
            charge = get_charge_attribute(atoms)
            spin_multiplicity = get_spin_multiplicity_attribute(atoms)

        if charge is not None:
            atoms.charge = charge  # type: ignore[attr-defined]
        if spin_multiplicity is not None:
            atoms.spin_multiplicity = spin_multiplicity  # type: ignore[attr-defined]

    # Strip the dummy atoms, if present
    if "X" in atoms.get_chemical_symbols():
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

    # Store Atoms object
    results["atoms"] = atoms

    return metadata | results | additional_fields
