"""Utility functions for dealing with defects."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires
from pymatgen.core.periodic_table import DummySpecies
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor

has_pmg_defects = bool(find_spec("pymatgen.analysis.defects"))
has_shakenbreak = bool(find_spec("shakenbreak"))
if has_pmg_defects:
    from pymatgen.analysis.defects.generators import VacancyGenerator
    from pymatgen.analysis.defects.thermo import DefectEntry
if has_shakenbreak:
    from shakenbreak.input import Distortions


if TYPE_CHECKING:
    from ase.atoms import Atoms
    from numpy.typing import NDArray
    from pymatgen.core.structure import Structure

    if has_pmg_defects:
        from pymatgen.analysis.defects.core import Defect
        from pymatgen.analysis.defects.generators import (
            AntiSiteGenerator,
            ChargeInterstitialGenerator,
            InterstitialGenerator,
            SubstitutionGenerator,
            VoronoiInterstitialGenerator,
        )


@requires(
    has_pmg_defects,
    "Missing pymatgen-analysis-defects. Please run pip install quacc[defects]",
)
@requires(has_shakenbreak, "Missing shakenbreak. Please run pip install quacc[defects]")
def make_defects_from_bulk(
    atoms: Atoms,
    defect_gen: (
        AntiSiteGenerator
        | ChargeInterstitialGenerator
        | InterstitialGenerator
        | SubstitutionGenerator
        | VacancyGenerator
        | VoronoiInterstitialGenerator
    ) = VacancyGenerator,
    defect_charge: int = 0,
    sc_mat: NDArray | None = None,
    min_atoms: int = 80,
    max_atoms: int = 240,
    min_length: float = 10.0,
    force_diagonal: bool = False,
    **defect_gen_kwargs,
) -> list[Atoms]:
    """
    Function to make defects from a bulk atoms object.

    Parameters
    ----------
    atoms
        bulk atoms
    defect_gen
        defect generator
    defect_charge
        charge state of defect
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
    **defect_gen_kwargs
        keyword arguments to pass to the pymatgen.analysis.defects.generators
        get_defects() method

    Returns
    -------
    list[Atoms]
        All generated defects
    """
    # Use pymatgen-analysis-defects and ShakeNBreak to generate defects
    struct = AseAtomsAdaptor.get_structure(atoms)

    # Make all the defects
    defects = defect_gen().get_defects(struct, **defect_gen_kwargs)
    final_defects = []
    for defect in defects:
        defect.user_charges = [defect_charge]

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
        defect_entry = get_defect_entry_from_defect(
            defect=defect,
            defect_supercell=defect_supercell,
            defect_charge=defect_charge,
        )

        # Instantiate class to apply rattle and bond distortion to all defects
        dist = Distortions([defect_entry])

        # Apply rattle and bond distortion to all defects
        defect_dict, distortion_metadata = dist.apply_distortions()
        defect_symbol = next(iter(distortion_metadata["defects"].keys()))
        distortion_dict = defect_dict[defect_symbol]["charges"][defect_charge][
            "structures"
        ]["distortions"]

        # Make atoms objects and store defect stats
        for distortions, defect_struct in distortion_dict.items():
            final_defect = defect_struct.to_ase_atoms()
            defect_stats = {
                "defect_symbol": defect_symbol,
                "defect_charge": defect_charge,
                "distortions": distortions,
                "bulk": atoms,
                "defect": defect,
            }
            final_defect.info["defect_stats"] = defect_stats
            final_defects.append(final_defect)
    return final_defects


@requires(
    has_pmg_defects,
    "Missing pymatgen-analysis-defects. Please run pip install quacc[defects]",
)
def get_defect_entry_from_defect(
    defect: Defect, defect_supercell: Structure, defect_charge: int
) -> DefectEntry:
    """
    Function to generate DefectEntry object from Defect object.

    Parameters
    ----------
    defect
        defect object
    defect_supercell
        defect supercell
    defect_charge
        charge state of defect

    Returns
    -------
    DefectEntry
        defect entry
    """
    # Find defect's fractional coordinates and remove it from supercell
    for site in defect_supercell:
        if site.species.elements[0].symbol == DummySpecies().symbol:
            dummy_site = site
            break
    defect_supercell.remove(dummy_site)

    computed_structure_entry = ComputedStructureEntry(
        structure=defect_supercell,
        energy=0.0,  # needs to be set, so set to 0.0
    )

    return DefectEntry(
        defect=defect,
        charge_state=defect_charge,
        sc_entry=computed_structure_entry,
        sc_defect_frac_coords=dummy_site.frac_coords,
    )
