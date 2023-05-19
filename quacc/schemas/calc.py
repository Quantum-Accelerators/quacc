"""
Schemas for storing ASE calculator data
"""
from __future__ import annotations

import os

from ase import Atoms
from ase.io import Trajectory
from atomate2.utils.path import get_uri

from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_


def summarize_run(
    atoms: Atoms,
    input_atoms: Atoms = None,
    prep_next_run: bool = True,
    additional_fields: dict = None,
) -> dict:
    """
    Get tabulated results from an Atoms object and calculator and store them in a database-friendly format.
    This is meant to be compatible with all calculator types.

    Parameters
    ----------
    atoms
        ASE Atoms following a calculation. A calculator must be attached.
    input_atoms
        Input ASE Atoms object to store.
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared for the next run
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    additional_fields
        Additional fields to add to the task document.

    Returns
    -------
    Dictionary representation of the task document with the following fields:
        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - database_version: str = Field(None, description="The database version for the built data.")
            - emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - pull_request: int = Field(None, description="The pull request number associated with this data build.")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - input_structure: Molecule | Structure = Field(None, title = "The Pymatgen Structure or Molecule object from the input Atoms object if input_atoms is not None.
            Inherits the Structure of Molecule metadata fields described below.
        - molecule: Molecule = Field(None, title = "The Pymatgen Molecule object from the Atoms object, if not periodic and store_pmg is True.")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - parameters: dict = Field(None, title = "The ase.calc.parameters used to run the calculation.")
        - structure: Structure = Field(None, title = "The Pymatgen Structure object from the Atoms object, if periodic and store_pmg is True.")
        - results: dict = Field(None, title = "The ase.calc.results from the calculation.")

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
    # Make sure there is a calculator with results
    if not atoms.calc:
        raise ValueError("ASE Atoms object has no attached calculator.")
    if not atoms.calc.results:
        raise ValueError("ASE Atoms object's calculator has no results.")

    if additional_fields is None:
        additional_fields = {}

    # Fetch all tabulated results from the attached calculator
    results = {"results": atoms.calc.results}

    # Get the calculator inputs
    uri = get_uri(os.getcwd())
    inputs = {
        "parameters": atoms.calc.parameters,
        "nid": uri.split(":")[0],
        "dir_name": ":".join(uri.split(":")[1:]),
    }
    if input_atoms:
        input_atoms_db = atoms_to_metadata(input_atoms)
        inputs["input_structure"] = input_atoms_db

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # Get tabulated properties of the structure itself
    atoms_db = atoms_to_metadata(atoms)

    # Create a dictionary of the inputs/outputs
    task_doc = {**atoms_db, **inputs, **results, **additional_fields}

    task_doc = dict(sorted(task_doc.items()))

    return task_doc


def summarize_opt_run(
    traj: Trajectory,
    parameters: dict,
    prep_next_run: bool = True,
    additional_fields: dict = None,
) -> dict:
    """
    Get tabulated results from an ASE Atoms trajectory and store them in a database-friendly format.
    This is meant to be compatible with all calculator types.

    Parameters
    ----------
    traj
        ASE Trajectory.
    parameters
        Dictionary of calculator parameters.
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared for the next run
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    additional_fields
        Additional fields to add to the task document.

    Returns
    -------
    Dict
        Dictionary of tabulated inputs/results

    """

    additional_fields = additional_fields or {}

    initial_atoms = traj[0]
    final_atoms = traj[-1]

    # Get results
    traj_results = {
        "trajectory_results": [atoms.calc.results for atoms in traj],
        "trajectory": [atoms_to_metadata(atoms) for atoms in traj],
    }
    results = {"results": final_atoms.calc.results}

    # Get the calculator inputs
    inputs = {
        "parameters": parameters,
        "dir_name": get_uri(os.getcwd()),
    }
    input_atoms_db = atoms_to_metadata(initial_atoms)
    inputs["input_structure"] = input_atoms_db

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        final_atoms = prep_next_run_(final_atoms)

    # Get tabulated properties of the structure itself
    atoms_db = atoms_to_metadata(final_atoms)

    # Create a dictionary of the inputs/outputs
    task_doc = {**atoms_db, **inputs, **results, **traj_results, **additional_fields}

    task_doc = dict(sorted(task_doc.items()))

    return task_doc
