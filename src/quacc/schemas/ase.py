"""Schemas for storing ASE-based data"""
from __future__ import annotations

import os
import warnings
from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.constraints import Filter
from ase.io import read
from ase.vibrations.data import VibrationsData

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.atoms import prep_next_run as prep_next_run_
from quacc.utils.db import results_to_db
from quacc.utils.dicts import clean_dict
from quacc.utils.files import get_uri

if TYPE_CHECKING:
    from ase import Atoms
    from ase.io import Trajectory
    from ase.optimize.optimize import Optimizer
    from ase.thermochemistry import IdealGasThermo
    from ase.vibrations import Vibrations
    from maggma.core import Store
    from tyipng import TypeVar

    RunSchema = TypeVar("RunSchema")
    OptSchema = TypeVar("OptSchema")
    VibSchema = TypeVar("VibSchema")
    ThermoSchema = TypeVar("ThermoSchema")


def summarize_run(
    atoms: Atoms,
    input_atoms: Atoms | None = None,
    charge_and_multiplicity: tuple[int, int] | None = None,
    prep_next_run: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> RunSchema:
    """
    Get tabulated results from an Atoms object and calculator and store them in a database-friendly format.
    This is meant to be compatible with all calculator types.

    Parameters
    ----------
    atoms
        ASE Atoms following a calculation. A calculator must be attached.
    input_atoms
        Input ASE Atoms object to store.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule metadata.
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared for the next run
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    remove_empties
        Whether to remove None values and empty lists/dicts from the task document.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None, `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:

        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - builder_meta.build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - builder_meta.database_version: str = Field(None, description="The database version for the built data.")
            - builder_meta.emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - builder_meta.pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - builder_meta.pull_request: int = Field(None, description="The pull request number associated with this data build.")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - input_structure: Molecule | Structure = Field(None, title = "The Pymatgen Structure or Molecule object from
         the input Atoms object if input_atoms is not None.")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - parameters: dict = Field(None, title = "the parameters used to run the calculation.")
        - results: dict = Field(None, title = "The results from the calculation.")

        For periodic structures, the task document also has the following fields:
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the material.")
        - composition: Composition = Field(None, description="Full composition for the material.")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition.")
        - density: float = Field(None, title="Density", description="Density in grams per cm^3.")
        - density_atomic: float = Field(None, title="Packing Density", description="The atomic packing density in
         atoms per cm^3.")
        - elements: List[Element] = Field(None, description="List of elements in the material.")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation
         of the formula.")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - nelements: int = Field(None, description="Number of elements.")
        - nsites: int = Field(None, description="Total number of sites in the structure.")
        - structure: Structure = Field(None, title = "The Pymatgen Structure object from the Atoms object, if
         periodic and store_pmg is True.")
        - symmetry: SymmetryData = Field(None, description="Symmetry data for this material.")
            - symmetry.crystal_system: CrystalSystem = Field(None, title="Crystal System", description="The crystal system for this lattice.")
            - symmetry.number: int = Field(None, title="Space Group Number", description="The spacegroup number for the lattice.")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice.")
            - symmetry.symbol: str = Field(None, title="Space Group Symbol", description="The spacegroup symbol for the lattice.")
            - symmetry.symprec: float = Field(None, title="Symmetry Finding Precision", description="The precision given to spglib to determine the symmetry of this lattice.")
        - volume: float = Field(None, title="Volume", description="Total volume for this structure in Angstroms^3.")

        For molecules that lack periodicity, the task document also has the following fields:

        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the
         molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular
         formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of
         the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - molecule: Molecule = Field(None, title = "The Pymatgen Molecule object from the Atoms object, if not
         periodic and store_pmg is True.")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons
         for the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - symmetry.eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - symmetry.linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - symmetry.matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")
            - symmetry.rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - symmetry.tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
    """
    # Make sure there is a calculator with results
    if not atoms.calc:
        msg = "ASE Atoms object has no attached calculator."
        raise ValueError(msg)
    if not atoms.calc.results:
        msg = "ASE Atoms object's calculator has no results."
        raise ValueError(msg)
    store = SETTINGS.PRIMARY_STORE if store is None else store

    additional_fields = additional_fields or {}

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
        input_atoms_db = atoms_to_metadata(
            input_atoms, charge_and_multiplicity=charge_and_multiplicity
        )
        inputs["input_atoms"] = input_atoms_db

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # Get tabulated properties of the structure itself
    atoms_db = atoms_to_metadata(atoms, charge_and_multiplicity=charge_and_multiplicity)

    # Create a dictionary of the inputs/outputs
    task_doc = clean_dict(
        atoms_db | inputs | results | additional_fields, remove_empties=remove_empties
    )

    if store:
        results_to_db(store, task_doc)

    return task_doc


def summarize_opt_run(
    dyn: Optimizer,
    trajectory: Trajectory | list[Atoms] = None,
    check_convergence: bool | None = None,
    charge_and_multiplicity: tuple[int, int] | None = None,
    prep_next_run: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> OptSchema:
    """
    Get tabulated results from an ASE Atoms trajectory and store them in a database-friendly format.
    This is meant to be compatible with all calculator types.

    Parameters
    ----------
    dyn
        ASE Optimizer object.
    trajectory
        ASE Trajectory object or list[Atoms] from reading a trajectory file.
        If None, the trajectory must be found in dyn.traj_atoms.
    check_convergence
        Whether to check the convergence of the calculation. Defaults to True in settings.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule metadata.
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared for the next run
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    remove_empties
        Whether to remove None values and empty lists/dicts from the task document.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None, `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:

        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - builder_meta.build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - builder_meta.database_version: str = Field(None, description="The database version for the built data.")
            - builder_meta.emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - builder_meta.pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - builder_meta.pull_request: int = Field(None, description="The pull request number associated with this data build.")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - input_structure: Molecule | Structure = Field(None, title = "The Pymatgen Structure or Molecule object from
         the input Atoms object if input_atoms is not None.")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - parameters: dict = Field(None, title = "the parameters used to run the calculation.")
        - opt_parameters: dict = Field(None, title = "the parameters used to run the optimization.")
        - results: dict = Field(None, title = "The results from the calculation.")
        - trajectory: List[Atoms] = Trajectory of Atoms objects
        - trajectory_results: List[dict] = List of ase.calc.results from the trajectory

        For periodic structures, the task document also has the following fields:

        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the material.")
        - composition: Composition = Field(None, description="Full composition for the material.")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition.")
        - density: float = Field(None, title="Density", description="Density in grams per cm^3.")
        - density_atomic: float = Field(None, title="Packing Density", description="The atomic packing density in
         atoms per cm^3.")
        - elements: List[Element] = Field(None, description="List of elements in the material.")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation
         of the formula.")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - nelements: int = Field(None, description="Number of elements.")
        - nsites: int = Field(None, description="Total number of sites in the structure.")
        - structure: Structure = Field(None, title = "The Pymatgen Structure object from the Atoms object, if periodic
         and store_pmg is True.")
        - symmetry: SymmetryData = Field(None, description="Symmetry data for this material.")
            - symmetry.crystal_system: CrystalSystem = Field(None, title="Crystal System", description="The crystal system for this lattice.")
            - symmetry.number: int = Field(None, title="Space Group Number", description="The spacegroup number for the lattice.")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice.")
            - symmetry.symbol: str = Field(None, title="Space Group Symbol", description="The spacegroup symbol for the lattice.")
            - symmetry.symprec: float = Field(None, title="Symmetry Finding Precision", description="The precision given to spglib to determine the symmetry of this lattice.")
        - volume: float = Field(None, title="Volume", description="Total volume for this structure in Angstroms^3.")

        For molecules that lack periodicity, the task document also has the following fields:

        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in
         the molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular
         formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation
         of the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - molecule: Molecule = Field(None, title = "The Pymatgen Molecule object from the Atoms object, if not
         periodic and store_pmg is True.")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons for
         the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - symmetry.eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - symmetry.linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - symmetry.matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")
            - symmetry.rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - symmetry.tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
    """

    check_convergence = (
        SETTINGS.CHECK_CONVERGENCE if check_convergence is None else check_convergence
    )
    additional_fields = additional_fields or {}
    opt_parameters = dyn.todict()
    if hasattr(dyn, "fmax"):
        opt_parameters = opt_parameters | {"fmax": dyn.fmax}
    store = SETTINGS.PRIMARY_STORE if store is None else store

    # Check convergence
    is_converged = dyn.converged()
    if check_convergence and not is_converged:
        msg = "Optimization did not converge."
        raise ValueError(msg)

    # Get trajectory
    if not trajectory:
        trajectory = (
            dyn.traj_atoms
            if hasattr(dyn, "traj_atoms")
            else read(dyn.trajectory.filename, index=":")
        )

    initial_atoms = trajectory[0]
    final_atoms = dyn.atoms.atoms if isinstance(dyn.atoms, Filter) else dyn.atoms

    # Get results
    traj_results = {
        "trajectory_results": [atoms.calc.results for atoms in trajectory],
        "trajectory": [
            atoms_to_metadata(atoms, charge_and_multiplicity=charge_and_multiplicity)
            for atoms in trajectory
        ],
    }

    results = {
        "results": final_atoms.calc.results
        | {"converged": is_converged, "nsteps": dyn.get_number_of_steps()}
    }

    # Get the calculator inputs
    uri = get_uri(os.getcwd())
    inputs = {
        "parameters": dyn.atoms.calc.parameters,
        "parameters_opt": opt_parameters,
        "nid": uri.split(":")[0],
        "dir_name": ":".join(uri.split(":")[1:]),
    }
    input_atoms_db = atoms_to_metadata(
        initial_atoms, charge_and_multiplicity=charge_and_multiplicity
    )
    inputs["input_structure"] = input_atoms_db

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        final_atoms = prep_next_run_(final_atoms)

    # Get tabulated properties of the structure itself
    atoms_db = atoms_to_metadata(
        final_atoms, charge_and_multiplicity=charge_and_multiplicity
    )

    # Create a dictionary of the inputs/outputs
    task_doc = clean_dict(
        atoms_db | inputs | results | traj_results | additional_fields,
        remove_empties=remove_empties,
    )

    if store:
        results_to_db(store, task_doc)

    return task_doc


def summarize_vib_run(
    vib: Vibrations,
    charge_and_multiplicity: tuple[int, int] | None = None,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> VibSchema:
    """
    Get tabulated results from an ASE Vibrations object and store them in a database-friendly format.

    Parameters
    ----------
    vib
        ASE Vibrations object.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule metadata.
    remove_empties
        Whether to remove None values and empty lists/dicts from the task document.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None, `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:

        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - builder_meta.build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - builder_meta.database_version: str = Field(None, description="The database version for the built data.")
            - builder_meta.emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - builder_meta.pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - builder_meta.pull_request: int = Field(None, description="The pull request number associated with this data build.")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - parameters: dict = Field(None, title = "the parameters used to run the calculation.")
        - vib_parameters: dict = Field(None, title = "the parameters used to run the vibrations.")
            - vib_parameters.delta: float = the Vibrations delta value
            - vib_parameters.direction: str = the Vibrations direction value
            - vib_parameters.method: str = the Vibrations method value
            - vib_parameters.ndof: int = the Vibrations ndof value
            - vib_parameters.nfree: int = the Vibrations nfree value
        - results: dict = Field(None, title = "The results from the calculation.")
            - results.imag_vib_freqs: List[float] = Imaginary vibrational frequencies in cm^-1
            - results.n_imag: int = Number of imaginary vibrational frequencies
            - results.vib_energies: List[float] = Vibrational energies in eV. 3N-5 or 3N-6 for molecules; 3N for solids.
            - results.vib_freqs: List[float] = Vibrational frequencies in cm^-1. 3N-5 or 3N-6 for molecules; 3N for solids.
            - results.vib_energies_raw: List[float] = Vibrational energies in eV of length 3N.
            - results.vib_freqs_raw: List[float] = Vibrational frequencies in cm^-1 of length 3N.

        For periodic structures, the task document also has the following fields:

        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the material.")
        - composition: Composition = Field(None, description="Full composition for the material.")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition.")
        - density: float = Field(None, title="Density", description="Density in grams per cm^3.")
        - density_atomic: float = Field(None, title="Packing Density", description="The atomic packing density in atoms
         per cm^3.")
        - elements: List[Element] = Field(None, description="List of elements in the material.")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the
         formula.")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - nelements: int = Field(None, description="Number of elements.")
        - nsites: int = Field(None, description="Total number of sites in the structure.")
        - structure: Structure = Field(None, title = "The Pymatgen Structure object from the Atoms object, if periodic
         and store_pmg is True.")
        - symmetry: SymmetryData = Field(None, description="Symmetry data for this material.")
            - symmetry.crystal_system: CrystalSystem = Field(None, title="Crystal System", description="The crystal system for this lattice.")
            - symmetry.number: int = Field(None, title="Space Group Number", description="The spacegroup number for the lattice.")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice.")
            - symmetry.symbol: str = Field(None, title="Space Group Symbol", description="The spacegroup symbol for the lattice.")
            - symmetry.symprec: float = Field(None, title="Symmetry Finding Precision", description="The precision given to spglib to determine the symmetry of this lattice.")
        - volume: float = Field(None, title="Volume", description="Total volume for this structure in Angstroms^3.")

        For molecules that lack periodicity, the task document also has the following fields:

        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the
         molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular
         formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of
         the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - molecule: Molecule = Field(None, title = "The Pymatgen Molecule object from the Atoms object, if not
         periodic and store_pmg is True.")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons for
         the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - symmetry.eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - symmetry.linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - symmetry.matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")
            - symmetry.rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - symmetry.tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
    """
    additional_fields = additional_fields or {}
    store = SETTINGS.PRIMARY_STORE if store is None else store

    vib_freqs_raw = vib.get_frequencies().tolist()
    vib_energies_raw = vib.get_energies().tolist()

    atoms = vib._atoms if isinstance(vib, VibrationsData) else vib.atoms

    # Convert imaginary modes to negative values for DB storage
    for i, f in enumerate(vib_freqs_raw):
        if np.imag(f) > 0:
            vib_freqs_raw[i] = -np.abs(f)
            vib_energies_raw[i] = -np.abs(vib_energies_raw[i])
        else:
            vib_freqs_raw[i] = np.abs(f)
            vib_energies_raw[i] = np.abs(vib_energies_raw[i])

    uri = get_uri(os.getcwd())
    inputs = {
        "parameters": None
        if isinstance(vib, VibrationsData)
        else atoms.calc.parameters,
        "parameters_vib": None
        if isinstance(vib, VibrationsData)
        else {
            "delta": vib.delta,
            "direction": vib.direction,
            "method": vib.method,
            "ndof": vib.ndof,
            "nfree": vib.nfree,
        },
        "nid": uri.split(":")[0],
        "dir_name": ":".join(uri.split(":")[1:]),
    }

    atoms_db = atoms_to_metadata(atoms, charge_and_multiplicity=charge_and_multiplicity)

    # Get the true vibrational modes
    natoms = len(atoms)
    if natoms == 1:
        vib_freqs = []
        vib_energies = []
    elif atoms.pbc.any():
        vib_freqs = vib_freqs_raw
        vib_energies = vib_energies_raw
    else:
        # Sort by absolute value
        vib_freqs_raw_sorted = vib_freqs_raw.copy()
        vib_energies_raw_sorted = vib_energies_raw.copy()
        vib_freqs_raw_sorted.sort(key=np.abs)
        vib_energies_raw_sorted.sort(key=np.abs)

        # Cut the 3N-5 or 3N-6 modes based on their absolute value
        n_modes = 3 * natoms - 5 if atoms_db["symmetry"]["linear"] else 3 * natoms - 6
        vib_freqs = vib_freqs_raw_sorted[-n_modes:]
        vib_energies = vib_energies_raw_sorted[-n_modes:]

    imag_vib_freqs = [f for f in vib_freqs if f < 0]

    results = {
        "results": {
            "imag_vib_freqs": imag_vib_freqs,
            "n_imag": len(imag_vib_freqs),
            "vib_energies": vib_energies,
            "vib_freqs": vib_freqs,
            "vib_energies_raw": vib_energies_raw,
            "vib_freqs_raw": vib_freqs_raw,
        }
    }

    task_doc = clean_dict(
        atoms_db | inputs | results | additional_fields, remove_empties=remove_empties
    )

    if store:
        results_to_db(store, task_doc)

    return task_doc


def summarize_thermo(
    igt: IdealGasThermo,
    temperature: float = 298.15,
    pressure: float = 1.0,
    charge_and_multiplicity: tuple[int, int] | None = None,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> ThermoSchema:
    """
    Get tabulated results from an ASE IdealGasThermo object and store them in a database-friendly format.

    Parameters
    ----------
    igt
        ASE IdealGasThermo object.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule metadata.
    remove_empties
        Whether to remove None values and empty lists/dicts from the task document.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None, `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:

        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - builder_meta.build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - builder_meta.database_version: str = Field(None, description="The database version for the built data.")
            - builder_meta.emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - builder_meta.pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - builder_meta.pull_request: int = Field(None, description="The pull request number associated with this data build.")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - thermo_parameters: dict = Field(None, title = "the parameters used to run the thermo calculation.")
            - thermo_parameters.temperature: float = Temperature in Kelvins
            - thermo_parameters.pressure: float = Pressure in bar
            - thermo_parameters.sigma: float = The rotational symmetry number of the molecule
            - thermo_parameters.spin_multiplicity: int = The spin multiplicity of the molecule
            - thermo_parameters.vib_freqs: List[float] = Vibrational frequencies in cm^-1 used for the thermo calculation
            - thermo_parameters.vib_energies: List[float] = Vibrational energies in eV used for the thermo calculation
            - thermo_parameters.n_imag: int = Number of imaginary vibrational frequencies ignored in the thermo calculation
        - results: dict = Field(None, title = "The results from the calculation.")
            - results.energy: float = The potential energy of the system in eV
            - results.enthalpy: float = The enthalpy of the system in eV
            - results.entropy: float = The entropy of the system in eV/K
            - results.gibbs_energy: float = The Gibbs free energy of the system in eV
            - results.zpe: float = The zero point vibrational energy of the system in eV

        The task document also has the following fields from the Molecule object:

        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the
         molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified
         representation of the composition")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular
         formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of
         the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the
         formula.")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - molecule: Molecule = Field(None, title = "The Pymatgen Molecule object from the Atoms object, if not
         periodic and store_pmg is True.")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons for
         the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - symmetry.eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - symmetry.linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - symmetry.matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")
            - symmetry.rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - symmetry.tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
    """

    additional_fields = additional_fields or {}
    store = SETTINGS.PRIMARY_STORE if store is None else store

    uri = get_uri(os.getcwd())
    spin_multiplicity = int(2 * igt.spin + 1)

    inputs = {
        "parameters_thermo": {
            "temperature": temperature,
            "pressure": pressure,
            "sigma": igt.sigma,
            "spin_multiplicity": spin_multiplicity,
            "vib_freqs": [e / units.invcm for e in igt.vib_energies],
            "vib_energies": igt.vib_energies.tolist(),
            "n_imag": igt.n_imag,
        },
        "nid": uri.split(":")[0],
        "dir_name": ":".join(uri.split(":")[1:]),
    }

    results = {
        "results": {
            "energy": igt.potentialenergy,
            "enthalpy": igt.get_enthalpy(temperature, verbose=True),
            "entropy": igt.get_entropy(temperature, pressure * 10**5, verbose=True),
            "gibbs_energy": igt.get_gibbs_energy(
                temperature, pressure * 10**5, verbose=True
            ),
            "zpe": igt.get_ZPE_correction(),
        }
    }

    if charge_and_multiplicity and spin_multiplicity != charge_and_multiplicity[1]:
        warnings.warn(
            "The IdealGasThermo spin multiplicity does not match the user-specified multiplicity.",
            UserWarning,
        )

    atoms_db = atoms_to_metadata(
        igt.atoms, charge_and_multiplicity=charge_and_multiplicity
    )

    task_doc = clean_dict(
        atoms_db | inputs | results | additional_fields, remove_empties=remove_empties
    )

    if store:
        results_to_db(store, task_doc)

    return task_doc
