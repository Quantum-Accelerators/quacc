"""Schemas for storing metadata about Atoms objects"""
from __future__ import annotations

from typing import Any, TypeVar

import numpy as np
from ase.atoms import Atom, Atoms
from emmet.core.structure import MoleculeMetadata, StructureMetadata
from monty.json import jsanitize
from pymatgen.io.ase import AseAtomsAdaptor

from quacc.util.atoms import copy_atoms
from quacc.util.dicts import clean_dict

AtomsSchema = TypeVar("AtomsSchema")


def atoms_to_metadata(
    atoms: Atoms,
    charge_and_multiplicity: tuple[int, int] | None = None,
    get_metadata: bool = True,
    strip_info: bool = False,
    store_pmg: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
) -> AtomsSchema:
    """
    Convert an ASE Atoms object to a dict suitable for storage in MongoDB.

    Parameters
    ----------
    atoms
        ASE Atoms object to store in {"atoms": atoms}
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule metadata.
    get_metadata
        Whether to store atoms metadata in the returned dict.
    strip_info
        Whether to strip the data from atoms.info in the returned {"atoms": atoms}.
        Note that this data will be stored in {"atoms_info": atoms.info} regardless
    store_pmg
        Whether to store the Pymatgen Structure/Molecule object in {"structure": Structure}
        or {"molecule": Molecule}, respectively.
    remove_empties
        Whether to remove None values and empty lists/dicts from the TaskDocument.
    additional_fields
        Additional fields to add to the document.

    Returns
    -------
    dict
        Dictionary representation of the atoms object with the following fields:
        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - database_version: str = Field(None, description="The database version for the built data.")
            - emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - pull_request: int = Field(None, description="The pull request number associated with this data build.")

        For periodic structures, the task document also has the following fields:
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the material.")
        - composition: Composition = Field(None, description="Full composition for the material.")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified representation of the composition.")
        - density: float = Field(None, title="Density", description="Density in grams per cm^3.")
        - density_atomic: float = Field(None, title="Packing Density", description="The atomic packing density in atoms per cm^3.")
        - elements: List[Element] = Field(None, description="List of elements in the material.")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the formula.")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the formula.")
        - nelements: int = Field(None, description="Number of elements.")
        - nsites: int = Field(None, description="Total number of sites in the structure.")
        - structure: Structure = Field(None, title = "The Pymatgen Structure object from the Atoms object, if periodic and store_pmg is True.")
        - symmetry: SymmetryData = Field(None, description="Symmetry data for this material.")
            - crystal_system: CrystalSystem = Field(None, title="Crystal System", description="The crystal system for this lattice.")
            - number: int = Field(None, title="Space Group Number", description="The spacegroup number for the lattice.")
            - point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice.")
            - symbol: str = Field(None, title="Space Group Symbol", description="The spacegroup symbol for the lattice.")
            - symprec: float = Field(None, title="Symmetry Finding Precision", description="The precision given to spglib to determine the symmetry of this lattice.")
        - volume: float = Field(None, title="Volume", description="Total volume for this structure in Angstroms^3.")

        For molecules that lack periodicity, the task document also has the following fields:
        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified representation of the composition")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the formula.")
        - molecule: Molecule = Field(None, title = "The Pymatgen Molecule object from the Atoms object, if not periodic and store_pmg is True.")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons for the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")
            - rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
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
            metadata = StructureMetadata().from_structure(struct).dict()
            if store_pmg:
                results["structure"] = struct
        else:
            mol = AseAtomsAdaptor().get_molecule(atoms, charge_spin_check=False)
            metadata = MoleculeMetadata().from_molecule(mol).dict()
            if store_pmg:
                results["molecule"] = mol
        metadata["builder_meta"]["build_date"] = str(
            metadata["builder_meta"]["build_date"]
        )
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
    atoms_doc = metadata | results | additional_fields

    return clean_dict(atoms_doc, remove_empties=remove_empties)


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
