"""Schemas for VASP"""
from __future__ import annotations

import os
import warnings
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, TypeVar

from emmet.core.tasks import TaskDoc
from maggma.core import Store
from pymatgen.command_line.bader_caller import bader_analysis_from_path
from pymatgen.command_line.chargemol_caller import ChargemolAnalysis

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_
from quacc.util.db import results_to_db
from quacc.util.dicts import clean_dict
from quacc.util.files import copy_decompress

if TYPE_CHECKING:
    from ase import Atoms


VaspSchema = TypeVar("VaspSchema")


def summarize_run(
    atoms: Atoms,
    dir_path: str | None = None,
    prep_next_run: bool = True,
    run_bader: bool | None = None,
    check_convergence: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> VaspSchema:
    """
    Get tabulated results from a VASP run and store them in a database-friendly format.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    dir_path
        Path to VASP outputs. A value of None specifies the current working directory
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared for the next run.
        This clears out any attached calculator and moves the final magmoms to the initial magmoms.
    run_bader
        Whether a Bader analysis should be performed. Will not run if bader executable is not in PATH even if
        bader is set to True.
        Defaults to VASP_BADER in settings.
    check_convergence
        Whether to throw an error if convergence is not reached.
    remove_empties
        Whether to remove None values and empty lists/dicts from the TaskDocument.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None, `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:

        - analysis: AnalysisDoc = Field(None, title="Calculation Analysis", description="Some analysis of calculation data after collection.")
            - analysis.delta_volume: float = Field(None, title="Volume Change", description="Volume change for the calculation.")
            - analysis.delta_volume_percent: float = Field(None, title="Volume Change Percent", description="Percent volume change for the calculation.")
            - analysis.max_force: float = Field(None, title="Max Force", description="Maximum force on any atom at the end of the calculation.")
            - analysis.warnings: List[str] = Field(None, title="Calculation Warnings",description="Warnings issued after analysis.")
            - analysis.errors: List[str] = Field(None, title="Calculation Errors", description="Errors issued after analysis.")
        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - builder_meta.build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - builder_meta.emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - builder_meta.pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the material.")
        - composition: Composition = Field(None, description="Full composition for the material.")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified representation of the composition.")
        - custodian: List[CustodianDoc] = Field(None, title="Calcs reversed data", description="Detailed custodian data for each VASP calculation contributing to the task document.")
            - custodian.corrections: List[Any] = Field(None, title="Custodian Corrections", description="List of custodian correction data for calculation.")
            - custodian.job: dict = Field(None, title="Cusotodian Job Data", description="Job data logged by custodian.",)
        - density: float = Field(None, title="Density", description="Density in grams per cm^3.")
        - density_atomic: float = Field(None, title="Packing Density", description="The atomic packing density in atoms per cm^3.")
        - dir_name: str = Field(None, description="The directory for this VASP task")
        - elements: List[Element] = Field(None, description="List of elements in the material.")
        - entry: ComputedEntry = Field(None, description="The ComputedEntry from the task doc")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the formula.")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the formula.")
        - input: InputDoc = Field(None, description="The input structure used to generate the current task document.")
            - input.incar: Dict[str, Any] = Field(None, description="INCAR parameters for the calculation")
            - input.is_hubbard: bool = Field(False, description="Is this a Hubbard +U calculation")
            - input.kpoints: Union[Dict[str, Any], Kpoints] = Field(None, description="KPOINTS for the calculation")
            - input.hubbards: Dict = Field(None, description="The hubbard parameters used")
            - input.lattice_rec: Lattice = Field(None, description="Reciprocal lattice of the structure")
            - input.nkpoints: int = Field(None, description="Total number of k-points")
            - input.potcar: List[str] = Field(None, description="POTCAR symbols in the calculation")
            - input.potcar_spec: List[PotcarSpec] = Field(None, description="Title and hash of POTCAR files used in the calculation")
            - input.potcar_type: List[str] = Field(None, description="List of POTCAR functional types.")
            - input.parameters: Dict = Field(None, description="Parameters from vasprun")
            - input.structure: Structure = Field(None, description="Input structure for the calculation")
        - nelements: int = Field(None, description="Number of elements.")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - nsites: int = Field(None, description="Total number of sites in the structure.")
        - orig_inputs: OrigInputs = Field(None, description="The exact set of input parameters used to generate the current task document.")
            - orig_inputs.incar: Union[Incar, Dict] = Field(None, description="Pymatgen object representing the INCAR file.")
            - orig_inputs.poscar: Poscar = Field(None, description="Pymatgen object representing the POSCAR file.")
            - orig_inputs.kpoints: Kpoints = Field(None, description="Pymatgen object representing the KPOINTS file.")
            - orig_inputs.potcar: Union[Potcar, VaspPotcar, List[Any]] = Field(None, description="Pymatgen object representing the POTCAR file.",)
        - output: OutputDoc = Field(None, description="The exact set of output parameters used to generate the current task document.")
            - output.bandgap: float = Field(None, description="The DFT bandgap for the last calculation")
            - output.density: float = Field(..., description="Density of in units of g/cc.")
            - output.direct_gap: float = the direct bandgap (eV)
            - output.dos_properties: DosProperties = Field(None, description="DOS properties for the material")
            - output.efermi: float = the fermi energy
            - output.energy: float = Field(..., description="Total Energy in units of eV.")
            - output.energy_per_atom: float = Field(None, description="The final DFT energy per atom for the last calculation")
            - output.forces: List[List[float]] = Field(None, description="The force on each atom in units of eV/A^2.")
            - output.ionic_steps: float = the number of ionic steps
            - output.is_gap_direct: bool = if the band gap is direct
            - output.mag_density: float = magnetization density
            - output.outcar: Outcar = Field(None, description="Pymatgen object representing the OUTCAR file.")
            - output.run_stats: Dict = Field(None, description="Runtime statistics from the calculation.")
            - output.stress: List[List[float]] = Field(None, description="The stress on the cell in units of kB.")
            - output.structure: Structure = Field(None, title="Output Structure", description="Output Structure from the VASP calculation.")
        - state: TaskState = Field(None, description="State of this calculation")
        - symmetry: SymmetryData = Field(None, description="Symmetry data for this material.")
            - symmetry.crystal_system: CrystalSystem = Field(None, title="Crystal System", description="The crystal system for this lattice.")
            - symmetry.number: int = Field(None, title="Space Group Number", description="The spacegroup number for the lattice.")
            - symmetry.point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice.")
            - symmetry.symbol: str = Field(None, title="Space Group Symbol", description="The spacegroup symbol for the lattice.")
            - symmetry.symprec: float = Field(None, title="Symmetry Finding Precision", description="The precision given to spglib to determine the symmetry of this lattice.")
            - symmetry.version
        - vasp_version: str: the version of VASP
        - volume: float = Field(None, title="Volume", description="Total volume for this structure in Angstroms^3.")

        If run_bader is True, the following fields are added:

        - bader
            - bader.atomic_volume: float = The atomic volume
            - bader.bader_charge: float = The net bader charge
            - bader.bader_spin: float = The net bader spin density
            - bader.bader_version: float = The bader version
            - bader.min_dist: List[float] = The bader min_dist parameter
            - bader.partial_charges: List[float] = The atom-projected bader partial charges
            - bader.spin_moments: List[float] = The atom-projected bader spin moments
    """

    additional_fields = additional_fields or {}
    run_bader = SETTINGS.VASP_BADER if run_bader is None else run_bader
    dir_path = dir_path or os.getcwd()
    store = SETTINGS.PRIMARY_STORE if store is None else store

    # Fetch all tabulated results from VASP outputs files
    # Fortunately, emmet already has a handy function for this
    results = TaskDoc.from_directory(dir_path).dict()
    uri = results["dir_name"]
    results["nid"] = uri.split(":")[0]
    results["dir_name"] = ":".join(uri.split(":")[1:])
    results["builder_meta"]["build_date"] = str(results["builder_meta"]["build_date"])

    # Check for calculation convergence
    if check_convergence and results["state"] != "successful":
        raise ValueError("VASP calculation did not converge. Will not store task data.")

    # Remove unnecessary fields
    for k in [
        "additional_json",
        "author",
        "calcs_reversed",
        "icsd_id",
        "last_updated",
        "structure",  # already in output
        "tags",
        "task_id",
        "task_label",
        "transformations",
        "vasp_objects",
    ]:
        results.pop(k, None)

    if "output" in results:
        results["output"].pop("elph_displaced_structures", None)
        results["output"].pop("frequency_dependent_dielectric", None)

    # Get Bader analysis
    if run_bader:
        try:
            bader_stats = bader_runner(dir_path)
        except Exception:
            bader_stats = None
            warnings.warn("Bader analysis could not be performed.", UserWarning)

        if bader_stats:
            results["bader"] = bader_stats

            # Attach bader charges/spins to structure object
            struct = results["output"]["structure"]
            struct.add_site_property("bader_charge", bader_stats["partial_charges"])
            if "spin_moments" in bader_stats:
                struct.add_site_property("bader_spin", bader_stats["spin_moments"])
            results["output"]["structure"] = struct

    # Prepares the Atoms object for the next run by moving the
    # final magmoms to initial, clearing the calculator state,
    # and assigning the resulting Atoms object a unique ID.
    if prep_next_run:
        atoms = prep_next_run_(atoms)

    # We use get_metadata=False and store_pmg=False because the TaskDocument already
    # makes the structure metadata for us
    atoms_db = atoms_to_metadata(atoms, get_metadata=False, store_pmg=False)

    # Make task document
    task_doc = clean_dict(
        results | atoms_db | additional_fields, remove_empties=remove_empties
    )

    # Store the results
    if store:
        results_to_db(store, task_doc)

    return task_doc


def bader_runner(path: str | None = None, scratch_dir: str | None = None) -> dict:
    """
    Runs a Bader partial charge and spin moment analysis using the VASP
    output files in the given path. This function requires that `bader`
    is located in your PATH environment variable. See
    http://theory.cm.utexas.edu/henkelman/code/bader for the bader code.

    Parameters
    ----------
    path
        The path where the VASP output files are located.
        Must include CHGCAR, AECCAR0, AECCAR2, and POTCAR files. These
        files can be gzip'd or not -- it doesn't matter.
        If None, the current working directory is used.
    scratch_dir
        The path where the Bader analysis will be run.
        Defaults to SETTINGS.SCRATCH_DIR.

    Returns
    -------
    dict
        Dictionary containing the Bader analysis summary:
            {
                "min_dist": List[float],
                "atomic_volume": List[float],
                "vacuum_charge": float,
                "vacuum_volume": float,
                "bader_version": float,
                "partial_charges": List[float],
                "spin_moments": List[float],
            }
    """
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    path = path or os.getcwd()

    # Make sure files are present
    relevant_files = ["AECCAR0", "AECCAR2", "CHGCAR", "POTCAR"]
    for f in relevant_files:
        if not os.path.exists(os.path.join(path, f)) and not os.path.exists(
            os.path.join(path, f"{f}.gz")
        ):
            raise FileNotFoundError(f"Could not find {f} in {path}.")

    # Run Bader analysis
    with TemporaryDirectory(dir=scratch_dir) as tmpdir:
        copy_decompress(relevant_files, tmpdir)
        bader_stats = bader_analysis_from_path(path)

    # Store the partial charge, which is much more useful than the
    # raw charge and is more intuitive than the charge transferred.
    # An atom with a positive partial charge is cationic, whereas
    # an atom with a negative partial charge is anionic.
    bader_stats["partial_charges"] = [-c for c in bader_stats["charge_transfer"]]

    # Some cleanup of the returned dictionary
    if "magmom" in bader_stats:
        bader_stats["spin_moments"] = bader_stats["magmom"]
    bader_stats.pop("charge", None)
    bader_stats.pop("charge_transfer", None)
    bader_stats.pop("reference_used", None)
    bader_stats.pop("magmom", None)

    return bader_stats


def chargemol_runner(
    path: str | None = None,
    atomic_densities_path: str | None = None,
    scratch_dir: str | None = None,
) -> dict:
    """
    Runs a Chargemol (i.e. DDEC6 + CM5) analysis using the VASP output files
    in the given path. This function requires that the chargemol executable,
    given by the name `Chargemol_09_26_2017_linux_parallel`,
    `Chargemol_09_26_2017_linux_serial`, or `chargemol` is in the system PATH
    environment variable. See https://sourceforge.net/projects/ddec/files for
    the Chargemol code.

    Parameters
    ----------
    path
        The path where the VASP output files are located.
        Must include CHGCAR, AECCAR0, AECCAR2, and POTCAR files. These
        files can be gzip'd or not -- it doesn't matter.
        If None, the current working directory is used.
    atomic_densities_path
        The path where the reference atomic densities are located for Chargemol.
        If None, we assume that this directory is defined in an environment variable
        named DDEC6_ATOMIC_DENSITIES_DIR.
        See the Chargemol documentation for more information.
    scratch_dir
        The path where the Chargemol analysis will be run.

    Returns
    -------
    dict
        Dictionary containing the Chargemol analysis summary:
            {
                "ddec": {
                            "partial_charges": List[float],
                            "spin_moments": List[float],
                            "dipoles": List[float],
                            "bond_order_sums": List[float],
                            "bond_order_dict": Dict
                        },
                "cm5": {
                            "partial_charges": List[float],
                        }
            }
    """
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    path = path or os.getcwd()

    # Make sure files are present
    relevant_files = ["AECCAR0", "AECCAR2", "CHGCAR", "POTCAR"]
    for f in relevant_files:
        if not os.path.exists(os.path.join(path, f)) and not os.path.exists(
            os.path.join(path, f"{f}.gz")
        ):
            raise FileNotFoundError(f"Could not find {f} in {path}.")

    # Check environment variable
    if atomic_densities_path is None and "DDEC6_ATOMIC_DENSITIES_DIR" not in os.environ:
        raise ValueError("DDEC6_ATOMIC_DENSITIES_DIR environment variable not defined.")

    # Run Chargemol analysis
    with TemporaryDirectory(dir=scratch_dir) as tmpdir:
        copy_decompress(relevant_files, tmpdir)
        chargemol_stats = ChargemolAnalysis(
            path=path,
            atomic_densities_path=atomic_densities_path,
        )

    # Some cleanup of the returned dictionary
    chargemol_stats.pop("rsquared_moments", None)
    chargemol_stats.pop("rcubed_moments", None)
    chargemol_stats.pop("rfourth_moments", None)

    return chargemol_stats
