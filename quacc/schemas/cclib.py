"""
Schemas for molecular DFT codes parsed by cclib
"""
from __future__ import annotations

import logging
import os
from typing import Type, TypeVar

import numpy as np
from ase import Atoms
from cclib.io import ccread
from cclib.method import (
    CSPA,
    DDEC6,
    LPA,
    MBO,
    MPA,
    Bader,
    Bickelhaupt,
    Density,
    Hirshfeld,
    volume,
)
from emmet.core.structure import MoleculeMetadata
from monty.json import jsanitize
from pymatgen.core import Molecule
from pymatgen.core.periodic_table import Element

from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_
from quacc.util.dicts import clean_dict
from quacc.util.files import find_recent_logfile, get_uri

__all__ = ["TaskDocument"]

_T = TypeVar("_T", bound="TaskDocument")


def summarize_run(
    atoms: Atoms,
    logfile_extensions: str | list[str],
    dir_path: str = None,
    pop_analysis: str | list[str] = None,
    check_convergence: bool = True,
    prep_next_run: bool = True,
    remove_empties: bool = False,
    additional_fields: dict = None,
) -> dict:
    """
    Get tabulated results from a molecular DFT run and store them in a database-friendly format.
    This is meant to be a general parser built on top of cclib.

    This document inherits from the following schemas: atomate2.common.schemas.cclib.TaskDocument,
    emmet.core.structure.MoleculeMetadata, emmet.core.base.EmmetBaseModel.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    logfile_extensions
        Possible extensions of the log file (e.g. ".log", ".out", ".txt", ".chk"). Note that
        only a partial match is needed. For instance, `.log` will match `.log.gz` and `.log.1.gz`.
        If multiple files with this extension are found, the one with the most recent change time
        will be used. For an exact match only, put in the full file name.
    dir_path
        The path to the folder containing the calculation outputs. A value of None specifies the
        current working directory.
    pop_analysis
        The name(s) of any cclib post-processing analysis to run. Note that for
        bader, ddec6, and hirshfeld, a cube file (.cube, .cub) must reside in dir_path.
        Supports: "cpsa", "mpa", "lpa", "bickelhaupt", "density", "mbo", "bader", "ddec6",
        "hirshfeld".
    check_convergence
         Whether to throw an error if geometry optimization convergence is not reached.
    prep_next_run
        Whether the Atoms object storeed in {"atoms": atoms} should be prepared for the next run.
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    remove_empties
        Whether to remove None values and empty lists/dicts from the TaskDocument.
    additional_fields
        Additional fields to add to the task document.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:
        - atoms: Atoms = Field(None, title = "The Atoms object obtained from atoms.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - attributes: Dict = Field(None, description="Computed properties and calculation outputs. See https://cclib.github.io/data_notes.html")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - database_version: str = Field(None, description="The database version for the built data.")
            - emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - pull_request: int = Field(None, description="The pull request number associated with this data build.")
        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified representation of the composition")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - energy: float = Field(None, description="Final total energy")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the formula.")
        - last_updated: str = Field(default_factory=datetime_str, description="Timestamp for this task document was last updated")
        - logfile: str = Field(None, description="Path to the log file used in the post-processing analysis")
        - metadata: Dict = Field(None, description="Calculation metadata, including input parameters and runtime statistics. See https://cclib.github.io/data_notes.html#metadata")
        - molecule: Molecule = Field(None, description="Final output molecule from the task")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons for the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - parameters: dict = Field(None, title = "The input parameters obtained from atoms.calc.parameters.")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")
            - point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
        - tags: List[str] = Field(None, description="Optional tags for this task document")
        - task_label: str = Field(None, description="A description of the task")
    """
    # Make sure there is a calculator with results
    if not atoms.calc:
        raise ValueError("ASE Atoms object has no attached calculator.")
    if not atoms.calc.results:
        raise ValueError("ASE Atoms object's calculator has no results.")

    additional_fields = additional_fields or {}
    dir_path = dir_path or os.getcwd()

    # Fortunately, there is already a cclib parser in Atomate2
    results = TaskDocument.from_logfile(
        dir_path, logfile_extensions, analysis=pop_analysis
    ).dict()
    uri = results["dir_name"]
    results["nid"] = uri.split(":")[0]
    results["dir_name"] = ":".join(uri.split(":")[1:])
    results["builder_meta"]["build_date"] = str(results["builder_meta"]["build_date"])
    results["logfile"] = results["logfile"].split(":")[-1]

    # Check convergence if requested
    if check_convergence and results["attributes"].get("optdone") is False:
        raise ValueError("Optimization not complete.")

    # Get the calculator inputs
    inputs = {"parameters": atoms.calc.parameters}

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # We use get_metadata=False and store_pmg=False because the TaskDocument already
    # makes the structure metadata for us
    atoms_db = atoms_to_metadata(atoms, get_metadata=False, store_pmg=False)

    # Create a dictionary of the inputs/outputs
    task_doc = {**atoms_db, **inputs, **results, **additional_fields}

    return clean_dict(task_doc, remove_empties=remove_empties)


class TaskDocument(MoleculeMetadata):
    """
    Definition of a cclib-generated task document.

    This can be used as a general task document for molecular DFT codes.
    For the list of supported packages, see https://cclib.github.io

    This is a mirror of atomate2.common.schemas.cclib.TaskDocument
    """

    molecule: Molecule = None
    energy: float = None
    dir_name: str = None
    logfile: str = None
    attributes: dict = None
    metadata: dict = None
    task_label: str = None
    tags: list[str] = None
    last_updated: str = None

    @classmethod
    def from_logfile(
        cls: Type[_T],
        dir_name: str,
        logfile_extensions: str | list[str],
        store_trajectory: bool = False,
        additional_fields: dict = None,
        analysis: str | list[str] = None,
        proatom_dir: str = None,
    ) -> _T:
        """
        Create a TaskDocument from a log file.

        For a full description of each field, see https://cclib.github.io/data.html.

        Parameters
        ----------
        dir_name
            The path to the folder containing the calculation outputs.
        logfile_extensions
            Possible extensions of the log file (e.g. ".log", ".out", ".txt", ".chk").
            Note that only a partial match is needed. For instance, `.log` will match
            `.log.gz` and `.log.1.gz`. If multiple files with this extension are found,
            the one with the most recent change time will be used. For an exact match
            only, put in the full file name.
        store_trajectory
            Whether to store the molecule objects along the course of the relaxation
            trajectory.
        additional_fields
            Dictionary of additional fields to add to TaskDocument.
        analysis
            The name(s) of any cclib post-processing analysis to run. Note that for
            bader, ddec6, and hirshfeld, a cube file (.cube, .cub) must be in dir_name.
            Supports: cpsa, mpa, lpa, bickelhaupt, density, mbo, bader, ddec6,
            hirshfeld.
        proatom_dir
            The path to the proatom directory if ddec6 or hirshfeld analysis are
            requested. See https://cclib.github.io/methods.html for details. If None,
            the PROATOM_DIR environment variable must point to the proatom directory.

        Returns
        -------
        TaskDocument
            A TaskDocument object summarizing the inputs/outputs of the log file.
        """

        # Find the most recent log file with the given extension in the
        # specified directory.
        logfile = find_recent_logfile(dir_name, logfile_extensions)
        if not logfile:
            raise FileNotFoundError(
                f"Could not find file with extension {logfile_extensions} in {dir_name}"
            )

        additional_fields = {} if additional_fields is None else additional_fields

        # Let's parse the log file with cclib
        cclib_obj = ccread(logfile)
        if not cclib_obj:
            raise ValueError(f"Could not parse {logfile}")

        # Fetch all the attributes (i.e. all input/outputs from cclib)
        attributes = jsanitize(cclib_obj.getattributes())

        # Store charge and multiplicity since we use it frequently
        charge = cclib_obj.charge
        mult = cclib_obj.mult

        # Let's move the metadata out of attributes for convenience and
        # store it separately
        attributes.pop("metadata")
        metadata = jsanitize(cclib_obj.metadata)

        # monty datetime bug workaround: github.com/materialsvirtuallab/monty/issues/275
        if metadata.get("wall_time", None):
            metadata["wall_time"] = [str(m) for m in metadata["wall_time"]]
        if metadata.get("cpu_time", None):
            metadata["cpu_time"] = [str(m) for m in metadata["cpu_time"]]

        # Get the final energy to store as its own key/value pair
        energy = (
            cclib_obj.scfenergies[-1] if cclib_obj.scfenergies is not None else None
        )

        # Now we construct the input molecule. Note that this is not necessarily
        # the same as the initial molecule from the relaxation because the
        # DFT package may have re-oriented the system. We only try to store
        # the input if it is XYZ-formatted though since the Molecule object
        # does not support internal coordinates or Gaussian Z-matrix.
        if (
            cclib_obj.metadata.get("coord_type", None) == "xyz"
            and cclib_obj.metadata.get("coords", None) is not None
        ):
            coords_obj = np.array(cclib_obj.metadata["coords"])
            input_species = [Element(e) for e in coords_obj[:, 0]]
            input_coords = coords_obj[:, 1:].tolist()
            input_molecule = Molecule(
                input_species,
                input_coords,
                charge=charge,
                spin_multiplicity=mult,
            )
            attributes["molecule_unoriented"] = input_molecule

        # These are duplicates of things made with MoleculeMetadata, so we
        # can just remove them here
        duplicates = ["atomnos", "atomcoords", "charge", "mult", "natom"]
        for duplicate in duplicates:
            attributes.pop(duplicate, None)

        # We will remove duplicates in the metadata too
        metadata_duplicates = ["coords", "coord_type"]
        for metadata_duplicate in metadata_duplicates:
            metadata.pop(metadata_duplicate, None)

        # Construct the Molecule object(s) from the trajectory
        species = [Element.from_Z(z) for z in cclib_obj.atomnos]
        coords = cclib_obj.atomcoords
        molecules = [
            Molecule(
                species,
                coord,
                charge=charge,
                spin_multiplicity=mult,
            )
            for coord in coords
        ]
        initial_molecule = molecules[0]
        final_molecule = molecules[-1]
        attributes["molecule_initial"] = initial_molecule
        if store_trajectory:
            attributes["trajectory"] = molecules

        # Store the HOMO/LUMO energies for convenience
        if cclib_obj.moenergies is not None and cclib_obj.homos is not None:
            homo_energies, lumo_energies, homo_lumo_gaps = _get_homos_lumos(
                cclib_obj.moenergies, cclib_obj.homos
            )
            attributes["homo_energies"] = homo_energies
            if lumo_energies:
                attributes["lumo_energies"] = lumo_energies
            if homo_lumo_gaps:
                attributes["homo_lumo_gaps"] = homo_lumo_gaps

                # The HOMO-LUMO gap for a spin-polarized system is ill-defined.
                # This is why we report both the alpha and beta channel gaps
                # above. Here, we report min(LUMO_alpha-HOMO_alpha,LUMO_beta-HOMO_beta)
                # in case the user wants to easily query by this too. For restricted
                # systems, this will always be the same as above.
                attributes["min_homo_lumo_gap"] = min(homo_lumo_gaps)

        # Calculate any properties
        if analysis:
            if type(analysis) is str:
                analysis = [analysis]
            analysis = [a.lower() for a in analysis]

            # Look for .cube or .cub files
            cubefile_path = find_recent_logfile(dir_name, [".cube", ".cub"])

            for analysis_name in analysis:
                if calc_attributes := cclib_calculate(
                    cclib_obj, analysis_name, cubefile_path, proatom_dir
                ):
                    attributes[analysis_name] = calc_attributes
                else:
                    attributes[analysis_name] = None

        doc = cls.from_molecule(
            final_molecule,
            energy=energy,
            dir_name=get_uri(dir_name),
            logfile=get_uri(logfile),
            attributes=attributes,
            metadata=metadata,
        )
        doc.molecule = final_molecule
        doc = doc.copy(update=additional_fields)
        return doc


def cclib_calculate(
    cclib_obj,
    method: str,
    cube_file: str,
    proatom_dir: str,
) -> dict | None:
    """
    Run a cclib population analysis.

    This is a mirror of atomate2.common.schemas.cclib.cclib_calculate

    Parameters
    ----------
    cclib_obj
        The cclib object to run the population analysis on.
    method
        The population analysis method to use.
    cube_file
        The path to the cube file to use for the population analysis.
        Needed only for Bader, DDEC6, and Hirshfeld
    proatom_dir
        The path to the proatom directory to use for the population analysis.
        Needed only for DDEC6 and Hirshfeld.

    Returns
    -------
    dict
        Dictionary of the calculated attributes.
    """

    method = method.lower()
    cube_methods = ["bader", "ddec6", "hirshfeld"]

    if method in cube_methods and not cube_file:
        raise FileNotFoundError(
            f"A cube file must be provided for {method}. Returning None."
        )
    if method in {"ddec6", "hirshfeld"} and not proatom_dir:
        if "PROATOM_DIR" not in os.environ:
            raise OSError("PROATOM_DIR environment variable not set. Returning None.")
        proatom_dir = os.path.expandvars(os.environ["PROATOM_DIR"])
    if proatom_dir and not os.path.exists(proatom_dir):
        raise FileNotFoundError(
            f"Protatom directory {proatom_dir} does not exist. Returning None."
        )

    if cube_file and method in cube_methods:
        vol = volume.read_from_cube(cube_file)

    if method == "bader":
        m = Bader(cclib_obj, vol)
    elif method == "bickelhaupt":
        m = Bickelhaupt(cclib_obj)
    elif method == "cpsa":
        m = CSPA(cclib_obj)
    elif method == "ddec6":
        m = DDEC6(cclib_obj, vol, str(proatom_dir))
    elif method == "density":
        m = Density(cclib_obj)
    elif method == "hirshfeld":
        m = Hirshfeld(cclib_obj, vol, str(proatom_dir))
    elif method == "lpa":
        m = LPA(cclib_obj)
    elif method == "mbo":
        m = MBO(cclib_obj)
    elif method == "mpa":
        m = MPA(cclib_obj)
    else:
        raise ValueError(f"{method} is not supported.")

    try:
        m.calculate()
    except AttributeError:
        return None

    # The list of available attributes after a calculation. This is hardcoded for now
    # until https://github.com/cclib/cclib/issues/1097 is resolved. Once it is, we can
    # delete this and just do `return calc_attributes.getattributes()`.
    avail_attributes = [
        "aoresults",
        "fragresults",
        "fragcharges",
        "density",
        "donations",
        "bdonations",
        "repulsions",
        "matches",
        "refcharges",
    ]
    return {
        attribute: getattr(m, attribute)
        for attribute in avail_attributes
        if hasattr(m, attribute)
    }


def _get_homos_lumos(
    moenergies: list[list[float]], homo_indices: list[int]
) -> tuple[list[float], list[float], list[float]]:
    """
    Calculate the HOMO, LUMO, and HOMO-LUMO gap energies in eV.

    This is a mirror of atomate2.common.schemas.cclib._get_homos_lumos

    Parameters
    ----------
    moenergies
        List of MO energies. For restricted calculations, list[list[float]] is
        length one. For unrestricted, it is length two.
    homo_indices
        Indices of the HOMOs.

    Returns
    -------
    homo_energies
        The HOMO energies (eV), split by alpha and beta
    lumo_energies
        The LUMO energies (eV), split by alpha and beta
    homo_lumo_gaps
        The HOMO-LUMO gaps (eV), calculated as LUMO_alpha-HOMO_alpha and
        LUMO_beta-HOMO_beta
    """
    homo_energies = [moenergies[i][h] for i, h in enumerate(homo_indices)]
    # Make sure that the HOMO+1 (i.e. LUMO) is in moenergies (sometimes virtual
    # orbitals aren't printed in the output)
    for i, h in enumerate(homo_indices):
        if len(moenergies[i]) < h + 2:
            return homo_energies, [], []
    lumo_energies = [moenergies[i][h + 1] for i, h in enumerate(homo_indices)]
    homo_lumo_gaps = [
        lumo_energies[i] - homo_energies[i] for i in range(len(homo_energies))
    ]
    return homo_energies, lumo_energies, homo_lumo_gaps
