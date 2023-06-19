"""
Schemas for VASP
"""
from __future__ import annotations

import os
import warnings

from ase import Atoms
from emmet.core.tasks import TaskDoc

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.util.atoms import prep_next_run as prep_next_run_
from quacc.util.dicts import clean_dict
from quacc.util.pop_analysis import bader_runner


def summarize_run(
    atoms: Atoms,
    dir_path: str | None = None,
    prep_next_run: bool = True,
    run_bader: bool = SETTINGS.VASP_BADER,
    check_convergence: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
) -> dict:
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
    check_convergence
        Whether to throw an error if convergence is not reached.
    remove_empties
        Whether to remove None values and empty lists/dicts from the TaskDocument.
    additional_fields
        Additional fields to add to the task document.

    Returns
    -------
    dict
        Dictionary representation of the task document with the following fields:

        The task document also inherits the following fields from emmet.core.TaskDoc:
        - analysis: AnalysisDoc = Field(None, title="Calculation Analysis", description="Some analysis of calculation data after collection.")
            - delta_volume: float = Field(None, title="Volume Change", description="Volume change for the calculation.")
            - delta_volume_percent: float = Field(None, title="Volume Change Percent", description="Percent volume change for the calculation.")
            - max_force: float = Field(None, title="Max Force", description="Maximum force on any atom at the end of the calculation.")
            - warnings: List[str] = Field(None, title="Calculation Warnings",description="Warnings issued after analysis.")
            - errors: List[str] = Field(None, title="Calculation Errors", description="Errors issued after analysis.")
        - atoms: Atoms = Field(None, title = "The Atoms object from the calculation result.")
        - atoms_info: dict = Field(None, title = "The Atoms object info obtained from atoms.info.")
        - builder_meta: EmmetMeta = Field(default_factory=EmmetMeta, description="Builder metadata."):
            - build_date: str = Field(default_factory=datetime.utcnow, description="The build date for this document.")
            - emmet_version: str = Field(__version__, description="The version of emmet this document was built with.")
            - pymatgen_version: str = Field(pmg_version, description="The version of pymatgen this document was built with.")
        - chemsys: str = Field(None, title="Chemical System", description="dash-delimited string of elements in the material.")
        - composition: Composition = Field(None, description="Full composition for the material.")
        - composition_reduced: Composition = Field(None, title="Reduced Composition", description="Simplified representation of the composition.")
        - custodian: List[CustodianDoc] = Field(None, title="Calcs reversed data", description="Detailed custodian data for each VASP calculation contributing to the task document.")
            - corrections: List[Any] = Field(None, title="Custodian Corrections", description="List of custodian correction data for calculation.")
            - job: dict = Field(None, title="Cusotodian Job Data", description="Job data logged by custodian.",)
        - density: float = Field(None, title="Density", description="Density in grams per cm^3.")
        - density_atomic: float = Field(None, title="Packing Density", description="The atomic packing density in atoms per cm^3.")
        - dir_name: str = Field(None, description="The directory for this VASP task")
        - elements: List[Element] = Field(None, description="List of elements in the material.")
        - entry: ComputedEntry = Field(None, description="The ComputedEntry from the task doc")
        - formula_anonymous: str = Field(None, title="Anonymous Formula", description="Anonymized representation of the formula.")
        - formula_pretty: str = Field(None, title="Pretty Formula", description="Cleaned representation of the formula.")
        - input: InputDoc = Field(None, description="The input structure used to generate the current task document.")
            - incar: Dict[str, Any] = Field(None, description="INCAR parameters for the calculation")
            - is_hubbard: bool = Field(False, description="Is this a Hubbard +U calculation")
            - kpoints: Union[Dict[str, Any], Kpoints] = Field(None, description="KPOINTS for the calculation")
            - hubbards: Dict = Field(None, description="The hubbard parameters used")
            - lattice_rec: Lattice = Field(None, description="Reciprocal lattice of the structure")
            - nkpoints: int = Field(None, description="Total number of k-points")
            - potcar: List[str] = Field(None, description="POTCAR symbols in the calculation")
            - potcar_spec: List[PotcarSpec] = Field(None, description="Title and hash of POTCAR files used in the calculation")
            - potcar_type: List[str] = Field(None, description="List of POTCAR functional types.")
            - parameters: Dict = Field(None, description="Parameters from vasprun")
            - structure: Structure = Field(None, description="Input structure for the calculation")
        - nelements: int = Field(None, description="Number of elements.")
        - nid: str = Field(None, title = "The node ID representing the machine where the calculation was run.")
        - nsites: int = Field(None, description="Total number of sites in the structure.")
        - orig_inputs: OrigInputs = Field(None, description="The exact set of input parameters used to generate the current task document.")
            - incar: Union[Incar, Dict] = Field(None, description="Pymatgen object representing the INCAR file.")
            - poscar: Poscar = Field(None, description="Pymatgen object representing the POSCAR file.")
            - kpoints: Kpoints = Field(None, description="Pymatgen object representing the KPOINTS file.")
            - potcar: Union[Potcar, VaspPotcar, List[Any]] = Field(None, description="Pymatgen object representing the POTCAR file.",)
        - output: OutputDoc = Field(None, description="The exact set of output parameters used to generate the current task document.")
            - bandgap: float = Field(None, description="The DFT bandgap for the last calculation")
            - density: float = Field(..., description="Density of in units of g/cc.")
            - direct_gap: float = the direct bandgap (eV)
            - dos_properties: DosProperties = Field(None, description="DOS properties for the material")
            - efermi: float = the fermi energy
            - energy: float = Field(..., description="Total Energy in units of eV.")
            - energy_per_atom: float = Field(None, description="The final DFT energy per atom for the last calculation")
            - forces: List[List[float]] = Field(None, description="The force on each atom in units of eV/A^2.")
            - ionic_steps: float = the number of ionic steps
            - is_gap_direct: bool = if the band gap is direct
            - mag_density: float = magnetization density
            - outcar: Outcar = Field(None, description="Pymatgen object representing the OUTCAR file.")
            - run_stats: Dict = Field(None, description="Runtime statistics from the calculation.")
            - stress: List[List[float]] = Field(None, description="The stress on the cell in units of kB.")
            - structure: Structure = Field(None, title="Output Structure", description="Output Structure from the VASP calculation.")
        - state: TaskState = Field(None, description="State of this calculation")
        - symmetry: SymmetryData = Field(None, description="Symmetry data for this material.")
            - crystal_system: CrystalSystem = Field(None, title="Crystal System", description="The crystal system for this lattice.")
            - number: int = Field(None, title="Space Group Number", description="The spacegroup number for the lattice.")
            - point_group: str = Field(None, title="Point Group Symbol", description="The point group for the lattice.")
            - symbol: str = Field(None, title="Space Group Symbol", description="The spacegroup symbol for the lattice.")
            - symprec: float = Field(None, title="Symmetry Finding Precision", description="The precision given to spglib to determine the symmetry of this lattice.")
            - version
        - vasp_version: str: the version of VASP
        - volume: float = Field(None, title="Volume", description="Total volume for this structure in Angstroms^3.")

        If run_bader is True, the following fields are added:
        - bader
            - atomic_volume: float = The atomic volume
            - bader_charge: float = The net bader charge
            - bader_spin: float = The net bader spin density
            - bader_version: float = The bader version
            - min_dist: List[float] = The bader min_dist parameter
            - partial_charges: List[float] = The atom-projected bader partial charges
            - spin_moments: List[float] = The atom-projected bader spin moments
    """

    additional_fields = additional_fields or {}
    dir_path = dir_path or os.getcwd()

    # Fetch all tabulated results from VASP outputs files
    # Fortunately, emmet already has a handy function for this
    results = TaskDoc.from_directory(dir_path).dict()
    uri = results["dir_name"]
    results["nid"] = uri.split(":")[0]
    results["dir_name"] = ":".join(uri.split(":")[1:])
    results["builder_meta"]["build_date"] = str(results["builder_meta"]["build_date"])

    # Check for calculation convergence
    if check_convergence and results["state"] != "successful":
        raise RuntimeError(
            "VASP calculation did not converge. Will not store task data."
        )

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

    task_doc = results | atoms_db | additional_fields

    return clean_dict(task_doc, remove_empties=remove_empties)
