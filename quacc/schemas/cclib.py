"""
Schemas for molecular DFT codes parsed by cclib
"""
from __future__ import annotations

import os

from ase import Atoms
from atomate2.common.schemas.cclib import TaskDocument

from quacc.schemas.atoms import _quacc_sanitize
from quacc.util.atoms import prep_next_run as prep_next_run_


def summarize_run(
    atoms: Atoms,
    logfile_extensions: str | list[str],
    dir_path: str = None,
    check_convergence: bool = True,
    transition_state: bool = False,
    prep_next_run: bool = True,
    additional_fields: dict = None,
) -> dict:
    """
    Get tabulated results from a molecular DFT run and store them in a database-friendly format.
    This is meant to be a general parser built on top of cclib.

    This document inherits from the following schemas:
    - atomate2.common.schemas.cclib.TaskDocument
        -> emmet.core.structure.MoleculeMetadata
            -> emmet.core.base.EmmetBaseModel

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
    check_convergence
         Whether to throw an error if convergence is not reached.
    transition_state
        Whether the calculation is a transition state (used for convergence check).
    prep_next_run
        Whether the Atoms object storeed in {"atoms": atoms} should be prepared for the next run.
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    additional_fields
        Additional fields to add to the task document.

    Returns
    -------
    Dictionary representation of the task document with the following fields:
        - atoms: Atoms = Field(None, title = "The Atoms object obtained from atoms.")
        - parameters: dict = Field(None, title = "The input parameters obtained from atoms.calc.parameters.")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")

        The task document also inherits the following fields from atomate2.common.schemas.cclib.TaskDocument:
        - molecule: Molecule = Field(None, description="Final output molecule from the task")
        - energy: float = Field(None, description="Final total energy")
        - dir_name: str = Field(None, description="Directory where the output is parsed")
        - logfile: str = Field(None, description="Path to the log file used in the post-processing analysis")
        - attributes: Dict = Field(None, description="Computed properties and calculation outputs. See https://cclib.github.io/data_notes.html")
        - metadata: Dict = Field(None, description="Calculation metadata, including input parameters and runtime statistics. See https://cclib.github.io/data_notes.html#metadata")
        - task_label: str = Field(None, description="A description of the task")
        - tags: List[str] = Field(None, description="Optional tags for this task document")
        - last_updated: str = Field(default_factory=datetime_str, description="Timestamp for this task document was last updated")
        - _schema: str = Field(__version__, description="Version of atomate2 used to create the document", alias="schema")

        ... as well as the following fields from emmet.core.structure.MoleculeMetadata:
        - charge: int = Field(None, description="Charge of the molecule")
        - spin_multiplicity: int = Field(None, description="Spin multiplicity of the molecule")
        - natoms: int = Field(None, description="Total number of atoms in the molecule")
        - elements: List[Element] = Field(None, description="List of elements in the molecule")
        - nelements: int = Field(None, title="Number of Elements")
        - nelectrons: int = Field(None, title="Number of electrons", description="The total number of electrons for the molecule")
        - composition: Composition = Field(None, description="Full composition for the molecule")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified representation of the composition")
        - formula_alphabetical: str = Field(None, title="Alphabetical Formula", description="Alphabetical molecular formula")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the formula.")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the formula")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the molecule")
        - symmetry: PointGroupData = Field(None, description="Symmetry data for this molecule")
            - point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice")
            - rotation_number: float = Field(None, title="Rotational Symmetry Number", description="Rotational symmetry number for the molecule")
            - linear: bool = Field(None, title="Molecule Linearity", description="Is the molecule linear?")
            - tolerance: float = Field(None, title="Point Group Analyzer Tolerance", description="Distance tolerance to consider sites as symmetrically equivalent.")
            - eigen_tolerance: float = Field(None, title="Interia Tensor Eigenvalue Tolerance", description="Tolerance to compare eigen values of the inertia tensor.")
            - matrix_tolerance: float = Field(None, title="Symmetry Operation Matrix Element Tolerance" description="Tolerance used to generate the full set of symmetry operations of the point group.")

        ... as well as the following fields from emmet.core.base.EmmetBaseModel:
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
            - pull_request: int = Field(None, description="The pull request number associated with this data build.")
            - database_version: str = Field(None, description="The database version for the built data.")
            - build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
    """
    # Make sure there is a calculator with results
    if not atoms.calc:
        raise ValueError("ASE Atoms object has no attached calculator.")
    if not atoms.calc.results:
        raise ValueError("ASE Atoms object's calculator has no results.")

    additional_fields = additional_fields or {}
    dir_path = dir_path or os.getcwd()

    # Fortunately, there is already a cclib parser in Atomate2
    results = TaskDocument.from_logfile(dir_path, logfile_extensions).dict()
    uri = results["dir_name"]
    results["nid"] = uri.split(":")[0]
    results["dir_name"] = ":".join(uri.split(":")[1:])
    results["builder_meta"]["build_date"] = str(results["builder_meta"]["build_date"])
    results["logfile"] = results["logfile"].split(":")[-1]

    # Check convergence if requested
    if check_convergence:
        # If it's an opt+freq job, we will just check convergence on the
        # frequency step. This is because sometimes the frequency job will
        # yield all positive modes, but the optimization tolerances may
        # not be entirely met. See https://gaussian.com/faq3.
        vibfreqs = results["attributes"].get("vibfreqs")
        if vibfreqs:
            n_imag = sum(vibfreq < 0 for vibfreq in vibfreqs)
            if n_imag >= 2:
                raise ValueError(f"Too many imaginary modes: {n_imag}")
            if n_imag == 1 and not transition_state:
                raise ValueError("One imaginary mode, but transition_state = False.")
            if n_imag == 0 and transition_state:
                raise ValueError("No imaginary modes, but transition_state = True.")
        elif results["attributes"].get("optdone") is False:
            raise ValueError("Optimization not complete.")

    # Get the calculator inputs
    inputs = {"parameters": atoms.calc.parameters}

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # Add atoms info
    results["atoms_info"] = _quacc_sanitize(atoms.info)

    # Create a dictionary of the inputs/outputs
    task_doc = {**inputs, **results, **additional_fields}

    # Sort dict
    task_doc = dict(sorted(task_doc.items()))

    return task_doc
