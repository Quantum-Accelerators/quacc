"""Schemas for molecular DFT codes parsed by cclib"""
from __future__ import annotations

import logging
import os
import warnings
from inspect import getmembers, isclass
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Literal, TypeVar

import cclib
from cclib.io import ccread
from emmet.core.structure import MoleculeMetadata
from maggma.core import Store
from monty.json import jsanitize
from pydantic import Field
from pymatgen.core import Molecule
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.atoms import prep_next_run as prep_next_run_
from quacc.utils.db import results_to_db
from quacc.utils.dicts import clean_dict
from quacc.utils.files import find_recent_logfile, get_uri

if TYPE_CHECKING:
    from ase import Atoms

    cclibSchema = TypeVar("cclibSchema")


def summarize_run(
    atoms: Atoms,
    logfile_extensions: str | list[str],
    dir_path: str | None = None,
    pop_analyses: list[
        Literal[
            "cpsa",
            "mpa",
            "lpa",
            "bickelhaupt",
            "density",
            "mbo",
            "bader",
            "ddec6",
            "hirshfeld",
        ]
    ]
    | None = None,
    check_convergence: bool | None = None,
    prep_next_run: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> cclibSchema:
    """
    Get tabulated results from a molecular DFT run and store them in a
    database-friendly format. This is meant to be a general parser built on top
    of cclib.

    This document inherits from the following schemas: cclibTaskDocument (see
    below), emmet.core.structure.MoleculeMetadata,
    emmet.core.base.EmmetBaseModel.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    logfile_extensions
        Possible extensions of the log file (e.g. ".log", ".out", ".txt",
        ".chk"). Note that only a partial match is needed. For instance, `.log`
        will match `.log.gz` and `.log.1.gz`. If multiple files with this
        extension are found, the one with the most recent change time will be
        used. For an exact match only, put in the full file name.
    dir_path
        The path to the folder containing the calculation outputs. A value of
        None specifies the current working directory.
    pop_analyses
        The name(s) of any cclib post-processing analysis to run. Note that for
        bader, ddec6, and hirshfeld, a cube file (.cube, .cub) must reside in
        dir_path. Supports: "cpsa", "mpa", "lpa", "bickelhaupt", "density",
        "mbo", "bader", "ddec6", "hirshfeld".
    check_convergence
         Whether to throw an error if geometry optimization convergence is not
         reached. Defaults to True in settings.
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared
        for the next run. This clears out any attached calculator and moves the
        final magmoms to the initial magmoms.
    remove_empties
        Whether to remove None values and empty lists/dicts from the
        TaskDocument.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following
        fields:

        - atoms: Atoms = Field(None, title = "The Atoms object obtained from
          atoms.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained
          from atoms.info.")
        - attributes: Dict = Field(None, description="Computed properties and
          calculation outputs. See https://cclib.github.io/data_notes.html")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta,
          description="Builder metadata."):
            - builder_meta.build_date: str =
              Field(default_factory=datetime.utcnow, description="The build date
              for this document.")
            - builder_meta.database_version: str = Field(None, description="The
              database version for the built data.")
            - builder_meta.emmet_version: str = Field(__version__,
              description="The version of emmet this document was built with.")
            - builder_meta.pymatgen_version: str = Field(pmg_version,
              description="The version of pymatgen this document was built
              with.")
            - builder_meta.pull_request: int = Field(None, description="The pull
              request number associated with this data build.")
        - charge: int = Field(None, description="Charge of the molecule")
        - chemsys: str = Field(None, title="Chemical System",
          description="dash-delimited string of elements in the molecule")
        - composition: Composition = Field(None, description="Full composition
          for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced
          Composition", description="Simplified representation of the
          composition")
        - dir_name: str = Field(None, description="Directory where the output is
          parsed")
        - elements: List[Element] = Field(None, description="List of elements in
          the molecule")
        - energy: float = Field(None, description="Final total energy")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula",
          description="Alphabetical molecular formula")
        - formula_anonymous: str = Field(None, title="Anonymous Formula",
          description="Anonymized representation of the formula")
        - formula_pretty: str = Field(None, title="Pretty Formula",
          description="Cleaned representation of the formula.")
        - logfile: str = Field(None, description="Path to the log file used in
          the post-processing analysis")
        - metadata: Dict = Field(None, description="Calculation metadata,
          including input parameters and runtime statistics. See
          https://cclib.github.io/data_notes.html#metadata")
        - molecule: Molecule = Field(None, description="Final output molecule
          from the task")
        - natoms: int = Field(None, description="Total number of atoms in the
          molecule")
        - nid: str = Field(None, title = "The node ID representing the machine
          where the calculation was run.")
        - nelectrons: int = Field(None, title="Number of electrons",
          description="The total number of electrons for the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - parameters: dict = Field(None, title = "The input parameters obtained
          from atoms.calc.parameters.")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of
          the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for
          this molecule")
            - symmetry.eigen_tolerance: float = Field(None, title="Interia
              Tensor Eigenvalue Tolerance", description="Tolerance to compare
              eigen values of the inertia tensor.")
            - symmetry.linear: bool = Field(None, title="Molecule Linearity",
              description="Is the molecule linear?")
            - symmetry.matrix_tolerance: float = Field(None, title="Symmetry
              Operation Matrix Element Tolerance" description="Tolerance used to
              generate the full set of symmetry operations of the point group.")
            - symmetry.point_group: str = Field(None, title="Point Group
              Symbol", description="The point group for the lattice")
            - symmetry.rotation_number: float = Field(None, title="Rotational
              Symmetry Number", description="Rotational symmetry number for the
              molecule")
            - symmetry.tolerance: float = Field(None, title="Point Group
              Analyzer Tolerance", description="Distance tolerance to consider
              sites as symmetrically equivalent.")
        - tags: List[str] = Field(None, description="Optional tags for this task
          document")
        - task_label: str = Field(None, description="A description of the task")
    """

    # Make sure there is a calculator with results
    if not atoms.calc:
        msg = "ASE Atoms object has no attached calculator."
        raise ValueError(msg)
    if not atoms.calc.results:
        msg = "ASE Atoms object's calculator has no results."
        raise ValueError(msg)
    store = SETTINGS.PRIMARY_STORE if store is None else store

    # Set defaults
    check_convergence = (
        SETTINGS.CHECK_CONVERGENCE if check_convergence is None else check_convergence
    )
    additional_fields = additional_fields or {}
    dir_path = dir_path or os.getcwd()

    # Fortunately, there is already a cclib parser in Atomate2
    results = _cclibTaskDocument.from_logfile(
        dir_path, logfile_extensions, store_trajectory=True, analysis=pop_analyses
    )
    uri = results["dir_name"]
    results["nid"] = uri.split(":")[0]
    results["dir_name"] = ":".join(uri.split(":")[1:])
    results["builder_meta"]["build_date"] = str(results["builder_meta"]["build_date"])
    results["logfile"] = results["logfile"].split(":")[-1]
    if results["attributes"].get("trajectory"):
        results["attributes"]["trajectory"] = [
            atoms_to_metadata(AseAtomsAdaptor().get_atoms(molecule))
            for molecule in results["attributes"]["trajectory"]
        ]

    # Check convergence if requested
    if check_convergence and results["attributes"].get("optdone") is False:
        msg = "Optimization not complete."
        raise ValueError(msg)

    # Get the calculator inputs
    inputs = {"parameters": atoms.calc.parameters}

    # Prepares the Atoms object for the next run by moving the final magmoms to
    # initial, clearing the calculator state, and assigning the resulting Atoms
    # object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # We use get_metadata=False and store_pmg=False because the TaskDocument
    # already makes the structure metadata for us
    atoms_db = atoms_to_metadata(atoms, get_metadata=False, store_pmg=False)

    # Create a dictionary of the inputs/outputs
    task_doc = clean_dict(
        atoms_db | inputs | results | additional_fields,
        remove_empties=remove_empties,
    )

    # Store the results
    if store:
        results_to_db(store, task_doc)

    return task_doc


class _cclibTaskDocument(MoleculeMetadata):
    """
    Definition of a cclib-generated task document.

    This can be used as a general task document for molecular DFT codes. For the
    list of supported packages, see https://cclib.github.io
    """

    molecule: Molecule = Field(None, description="Final output molecule from the task")
    energy: float = Field(None, description="Final total energy")
    dir_name: str = Field(None, description="Directory where the output is parsed")
    logfile: str = Field(
        None, description="Path to the log file used in the post-processing analysis"
    )
    attributes: Dict = Field(
        None, description="Computed properties and calculation outputs"
    )
    metadata: Dict = Field(
        None,
        description="Calculation metadata, including input parameters and runtime "
        "statistics",
    )
    task_label: str = Field(None, description="A description of the task")
    tags: List[str] = Field(None, description="Optional tags for this task document")

    @classmethod
    def from_logfile(
        cls,
        dir_name: str | Path,
        logfile_extensions: str | list[str],
        store_trajectory: bool = False,
        additional_fields: dict | None = None,
        analysis: str | list[str] | None = None,
        proatom_dir: Path | str | None = None,
    ) -> dict:
        """
        Create a TaskDocument from a log file.

        For a full description of each field, see
        https://cclib.github.io/data.html.

        This is a mirror of atomate2.common.schemas.cclib.TaskDocument

        Parameters
        ----------
        dir_name
            The path to the folder containing the calculation outputs.
        logfile_extensions
            Possible extensions of the log file (e.g. ".log", ".out", ".txt",
            ".chk"). Note that only a partial match is needed. For instance,
            `.log` will match `.log.gz` and `.log.1.gz`. If multiple files with
            this extension are found, the one with the most recent change time
            will be used. For an exact match only, put in the full file name.
        store_trajectory
            Whether to store the molecule objects along the course of the
            relaxation trajectory.
        additional_fields
            Dictionary of additional fields to add to TaskDocument.
        analysis
            The name(s) of any cclib post-processing analysis to run. Note that
            for bader, ddec6, and hirshfeld, a cube file (.cube, .cub) must be
            in dir_name. Supports: cpsa, mpa, lpa, bickelhaupt, density, mbo,
            bader, ddec6, hirshfeld.
        proatom_dir
            The path to the proatom directory if ddec6 or hirshfeld analysis are
            requested. See https://cclib.github.io/methods.html for details. If
            None, the PROATOM_DIR environment variable must point to the proatom
            directory.

        Returns
        -------
        dict
            A TaskDocument dictinoary summarizing the inputs/outputs of the log
            file.
        """

        # Find the most recent log file with the given extension in the
        # specified directory.
        logfile = find_recent_logfile(dir_name, logfile_extensions)
        if not logfile:
            msg = (
                f"Could not find file with extension {logfile_extensions} in {dir_name}"
            )
            raise FileNotFoundError(msg)

        additional_fields = {} if additional_fields is None else additional_fields

        # Let's parse the log file with cclib
        cclib_obj = ccread(logfile, logging.ERROR)
        if not cclib_obj:
            msg = f"Could not parse {logfile}"
            raise ValueError(msg)

        # Fetch all the attributes (i.e. all input/outputs from cclib)
        attributes = jsanitize(cclib_obj.getattributes())

        # Store charge and multiplicity since we use it frequently
        charge = cclib_obj.charge
        mult = cclib_obj.mult

        # Let's move the metadata out of attributes for convenience and store it
        # separately
        attributes.pop("metadata")
        metadata = jsanitize(cclib_obj.metadata)

        # monty datetime bug workaround:
        # github.com/materialsvirtuallab/monty/issues/275
        if wall_time := metadata.get("wall_time"):
            metadata["wall_time"] = [*map(str, wall_time)]
        if cpu_time := metadata.get("cpu_time"):
            metadata["cpu_time"] = [*map(str, cpu_time)]

        # Get the final energy to store as its own key/value pair
        energy = (
            cclib_obj.scfenergies[-1] if cclib_obj.scfenergies is not None else None
        )

        # Now we construct the input molecule. Note that this is not necessarily
        # the same as the initial molecule from the relaxation because the DFT
        # package may have re-oriented the system. We only try to store the
        # input if it is XYZ-formatted though since the Molecule object does not
        # support internal coordinates or Gaussian Z-matrix.
        if (
            cclib_obj.metadata.get("coord_type") == "xyz"
            and cclib_obj.metadata.get("coords") is not None
        ):
            coords_obj = cclib_obj.metadata["coords"]
            input_species = [Element(row[0]) for row in coords_obj]
            input_coords = [row[1:] for row in coords_obj]
            input_molecule = Molecule(
                input_species,
                input_coords,
                charge=charge,
                spin_multiplicity=mult,
            )
            attributes["molecule_unoriented"] = input_molecule

        # These are duplicates of things made with MoleculeMetadata, so we can
        # just remove them here
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
                # above. Here, we report
                # min(LUMO_alpha-HOMO_alpha,LUMO_beta-HOMO_beta) in case the
                # user wants to easily query by this too. For restricted
                # systems, this will always be the same as above.
                attributes["min_homo_lumo_gap"] = min(homo_lumo_gaps)

        # Calculate any properties
        if analysis:
            if isinstance(analysis, str):
                analysis = [analysis]
            analysis = [a.lower() for a in analysis]

            # Look for .cube or .cub files
            cubefile_path = find_recent_logfile(dir_name, [".cube", ".cub"])

            for analysis_name in analysis:
                if calc_attributes := _cclib_calculate(
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
        return doc.copy(update=additional_fields).dict()


def _cclib_calculate(
    cclib_obj,
    method: str,
    cube_file: Path | str | None = None,
    proatom_dir: Path | str | None = None,
) -> dict | None:
    """
    Run a cclib population analysis.

    Parameters
    ----------
    cclib_obj
        The cclib object to run the population analysis on.
    method
        The population analysis method to use.
    cube_file
        The path to the cube file to use for the population analysis. Needed
        only for Bader, DDEC6, and Hirshfeld
    proatom_dir
        The path to the proatom directory to use for the population analysis.
        Needed only for DDEC6 and Hirshfeld.
    """

    method = method.lower()
    cube_methods = ["bader", "ddec6", "hirshfeld"]
    proatom_methods = ["ddec6", "hirshfeld"]

    if method in cube_methods:
        if not cube_file:
            msg = f"A cube file must be provided for {method}."
            raise ValueError(msg)
        if not os.path.exists(cube_file):
            msg = f"Cube file {cube_file} does not exist."
            raise FileNotFoundError(msg)
    if method in proatom_methods:
        if not proatom_dir:
            if os.getenv("PROATOM_DIR") is None:
                msg = "PROATOM_DIR environment variable or proatom_dir kwarg needs to be set."
                raise ValueError(msg)
            proatom_dir = os.path.expandvars(os.environ["PROATOM_DIR"])
        if not os.path.exists(proatom_dir):
            msg = f"Protatom directory {proatom_dir} does not exist. Returning None."
            raise FileNotFoundError(msg)
    cclib_methods = getmembers(cclib.method, isclass)
    method_class = next(
        (
            cclib_method[1]
            for cclib_method in cclib_methods
            if cclib_method[0].lower() == method
        ),
        None,
    )
    if method_class is None:
        msg = f"{method} is not a valid cclib population analysis method."
        raise ValueError(msg)

    if method in cube_methods:
        vol = cclib.method.volume.read_from_cube(str(cube_file))
        if method in proatom_methods:
            m = method_class(cclib_obj, vol, str(proatom_dir))
        else:
            m = method_class(cclib_obj, vol)
    else:
        m = method_class(cclib_obj)

    try:
        m.calculate()
    except Exception as e:
        warnings.warn(f"Could not calculate {method}: {e}", UserWarning)
        return None

    # The list of available attributes after a calculation. This is hardcoded
    # for now until https://github.com/cclib/cclib/issues/1097 is resolved. Once
    # it is, we can delete this and just do `return
    # calc_attributes.getattributes()`.
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
) -> tuple[list[float], list[float], list[float]] | tuple[list[float], None, None]:
    """
    Calculate the HOMO, LUMO, and HOMO-LUMO gap energies in eV.

    Parameters
    ----------
    moenergies
        List of MO energies. For restricted calculations, List[List[float]] is
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
            return homo_energies, None, None
    lumo_energies = [moenergies[i][h + 1] for i, h in enumerate(homo_indices)]
    homo_lumo_gaps = [
        lumo_energies[i] - homo_energies[i] for i in range(len(homo_energies))
    ]
    return homo_energies, lumo_energies, homo_lumo_gaps
