"""
Utility functions for dealing with defects
"""
from __future__ import annotations

import numpy as np
from ase.atoms import Atoms
from pymatgen.analysis.defects.core import Defect
from pymatgen.analysis.defects.generators import (
    AntiSiteGenerator,
    ChargeInterstitialGenerator,
    InterstitialGenerator,
    SubstitutionGenerator,
    VacancyGenerator,
    VoronoiInterstitialGenerator,
)
from pymatgen.analysis.defects.thermo import DefectEntry
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecies
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor
from shakenbreak.input import Distortions

# NOTES:
# - Anytime an Atoms object is converted to a pmg structure, make sure
# to reattach any .info flags to the Atoms object, e.g. via `new_atoms.info = atoms.info.copy()``.
# Note that atoms.info is mutable, so copy it!
# - All major functions should take in Atoms by default and return Atoms
# by default. Pymatgen structures can be returned with an optional kwarg.
# - If you modify the properties of an input Atoms object in any way, make sure to do so
# on a copy because Atoms objects are mutable.


def make_defects_from_bulk(
    atoms: Atoms,
    defectgen: (
        AntiSiteGenerator
        | ChargeInterstitialGenerator
        | InterstitialGenerator
        | SubstitutionGenerator
        | VacancyGenerator
        | VoronoiInterstitialGenerator
    ) = VacancyGenerator,
    charge_state: int = 0,
    sc_mat: np.ndarray | None = None,
    min_atoms: int = 80,
    max_atoms: int = 240,
    min_length: float = 10.0,
    force_diagonal: bool = False,
    **defectgen_kwargs,
) -> list[Atoms]:
    """
    Function to make defects from a bulk atoms object.

    Parameters
    ----------
    atoms
        bulk atoms
    defectgen
        defect generator
    sc_mat
        supercell matrix
    min_atoms
        minimum number of atoms in supercell
    max_atoms
        maximum number of atoms in supercell
    min_length
        minimum length of supercell
    force_diagonal
        force supercell to be diagonal
    charge_state
        charge state of defect
    **defectgen_kwargs
        keyword arguments to pass to the pymatgen.analysis.defects.generators get_defects() method

    Returns
    -------
    list[Atoms]
        All generated defects
    """

    # Use pymatgen-analysis-defects and ShakeNBreak to generate defects
    struct = AseAtomsAdaptor.get_structure(atoms)
    atoms_info = atoms.info.copy()

    # Make all the defects
    defects = defectgen().get_defects(struct, **defectgen_kwargs)
    final_defects = []
    for defect in defects:
        defect.user_charges = [charge_state]

        # Generate the supercell for a defect
        defect_supercell = defect.get_supercell_structure(
            sc_mat=sc_mat,
            dummy_species=DummySpecies(),
            min_atoms=min_atoms,
            max_atoms=max_atoms,
            min_length=min_length,
            force_diagonal=force_diagonal,
        )

        # Generate DefectEntry object from Defect object
        defect_entry = _get_defect_entry_from_defect(
            defect=defect, defect_supercell=defect_supercell, charge_state=charge_state
        )

        # Instantiate class to apply rattle and bond distortion to all defects
        Dist = Distortions(defects=[defect_entry])

        # Apply rattle and bond distortion to all defects
        defect_dict, distortion_metadata = Dist.apply_distortions()
        defect_symbol = list(distortion_metadata["defects"].keys())[0]
        distortion_dict = defect_dict[defect_symbol]["charges"][charge_state][
            "structures"
        ]["distortions"]

        # Make atoms objects and store defect stats
        for distortions, defect_struct in distortion_dict.items():
            final_defect = AseAtomsAdaptor.get_atoms(defect_struct)
            final_defect.info = atoms_info.copy()
            defect_stats = {
                "defect_symbol": defect_symbol,
                "charge_state": charge_state,
                "distortions": distortions,
                "bulk": atoms,
                "defect": defect,
            }
            final_defect.info["defect_stats"] = defect_stats
            final_defects.append(final_defect)
    return final_defects


def _get_defect_entry_from_defect(
    defect: Defect,
    defect_supercell: Structure,
    charge_state: int,
) -> DefectEntry:
    """
    Function to generate DefectEntry object from Defect object

    Parameters
    ----------
    defect
        defect object
    defect_supercell
        defect supercell
    charge_state
        charge state of defect

    Returns
    -------
    DefectEntry
        defect entry
    """
    # Find defect's fractional coordinates and remove it from supercell
    dummy_site = [
        site
        for site in defect_supercell
        if site.species.elements[0].symbol == DummySpecies().symbol
    ][0]
    sc_defect_frac_coords = dummy_site.frac_coords
    defect_supercell.remove(dummy_site)

    computed_structure_entry = ComputedStructureEntry(
        structure=defect_supercell,
        energy=0.0,  # needs to be set, so set to 0.0
    )

    return DefectEntry(
        defect=defect,
        charge_state=charge_state,
        sc_entry=computed_structure_entry,
        sc_defect_frac_coords=sc_defect_frac_coords,
    )
