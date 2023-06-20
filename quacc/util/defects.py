"""
Utility functions for dealing with defects
"""
from __future__ import annotations

from ase.atoms import Atoms
from pymatgen.analysis.defects.core import Defect
from pymatgen.analysis.defects.generators import AntiSiteGenerator, ChargeInterstitialGenerator, \
    InterstitialGenerator, SubstitutionGenerator, VacancyGenerator, VoronoiInterstitialGenerator
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


def get_defect_entry_from_defect(
        defect: Defect,
        defect_supercell: Structure,
        charge_state: int,
        dummy_species: DummySpecies,
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
    dummy_species
        Dummy species (used to keep track of the defect coords in the supercell)

    Returns
    -------
    DefectEntry
        defect entry
    """
    # Find defect's fractional coordinates and remove it from supercell
    dummy_site = [
        site for site in defect_supercell
        if site.species.elements[0].symbol == dummy_species.symbol
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


def make_defects_from_bulk(
        atoms: Atoms,
        def_gen: (AntiSiteGenerator | ChargeInterstitialGenerator | InterstitialGenerator | SubstitutionGenerator
                  | VacancyGenerator | VoronoiInterstitialGenerator),
        sc_mat: 'np.ndarray | None' = None,
        dummy_species: DummySpecies = DummySpecies("X"),
        min_atoms: int = 80,
        max_atoms: int = 240,
        min_length: float = 10.0,
        force_diagonal: bool = False,
        charge_state: int = 0,
        **defectgen_kwargs,  # TODO: what are these?
) -> list[Atoms]:
    """
    Function to make defects from a bulk atoms object.

    Parameters
    ----------
    atoms
        bulk atoms
    def_gen
        defect generator
    sc_mat
        supercell matrix
    dummy_species
        A special specie for representing non-traditional elements or species. For example, representation of vacancies
        (charged or otherwise), or special sites, etc.
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
    **defectgen_kwargs: keyword arguments to pass to  # TODO: what are these?

    Returns
    -------
    list[Atoms]
        All generated defects
    """

    # Use pymatgen-analysis-defects and ShakeNBreak to generate defects
    struct = AseAtomsAdaptor.get_structure(atoms)
    defects = def_gen().get_defects(struct)

    # TODO: select defects
    # Make all the defects
    defect_structs = []  # TODO: unnecessary?
    # for defect in defects:
    for defect in [defects[1]]:  # TODO: remove
        defect.user_charges = [charge_state]

        # Generate the supercell for a defect
        defect_supercell = defect.get_supercell_structure(
            sc_mat=sc_mat,
            dummy_species=dummy_species,
            min_atoms=min_atoms,
            max_atoms=max_atoms,
            min_length=min_length,
            force_diagonal=force_diagonal,
        )

        # Generate DefectEntry object from Defect object
        defect_entry = get_defect_entry_from_defect(
            defect=defect,
            defect_supercell=defect_supercell,
            charge_state=charge_state,
            dummy_species=dummy_species,
        )

        # Instantiate class to apply rattle and bond distortion to all defects
        Dist = Distortions(defects=[defect_entry])

        # Apply rattle and bond distortion to all defects
        defect_dict, distortion_metadata = Dist.apply_distortions()
        defect_symbol = list(distortion_metadata["defects"].keys())[0]
        charge_state = list(distortion_metadata["defects"][defect_symbol]["charges"].keys())[0]
        defect_structs = defect_dict[defect_symbol]["charges"][charge_state]["structures"]["distortions"].values()

        # Make atoms objects
        # TODO: add distortion metadata to atoms.info?
        final_defects = [AseAtomsAdaptor.get_atoms(defect_struct) for defect_struct in defect_structs]
        break  # TODO: remove
    return final_defects
