from ase.atoms import Atoms
import ase.io.jsonio as jsonio
from ase.calculators.singlepoint import SinglePointDFTCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import generate_all_slabs, Slab
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
import numpy as np

# properties supported by SinglePointDFTCalculator
ALL_PROPERTIES = [
    "energy",
    "forces",
    "stress",
    "stresses",
    "dipole",
    "charges",
    "magmom",
    "magmoms",
    "free_energy",
    "energies",
]


def serialize(atoms):
    """
    A smarter version of ASE's encode() function for Atoms. Stores the calculator
    results in atoms.info["results"] for later retrieval. This makes it
    possible to do things like atoms.get_magnetic_moment() even after
    an encode/decode cycle.

    Args:
        atoms (ase.Atoms): Atoms object

    Returns:
        atoms (ase.Atoms): Atoms object with calculator results attached in atoms.info["results"]
    """
    if getattr(atoms, "calc", None) and getattr(atoms.calc, "results", None):
        atoms.info["results"] = atoms.calc.results
    atoms = jsonio.encode(atoms)
    return atoms


def deserialize(atoms):
    """
    A smarter version of ASE's decode() function for Atoms. Retrieves the calculator
    results from atoms.info["results"] to re-instantiate. This makes it
    possible to do things like atoms.get_magnetic_moment() even after
    an encode/decode cycle.

    Args:
        atoms (ase.Atoms): Atoms object

    Returns:
        atoms (ase.Atoms): Atoms object with calculator results attached in atoms.info["results"]
    """
    atoms = jsonio.decode(atoms)

    supported_results = {}

    # Store results in a SinglePointDFTCalculator to make sure they stay between
    # serialize/deserialize cycles
    if atoms.info.get("results", None):
        for prop, result in atoms.info["results"].items():
            if prop in ALL_PROPERTIES:
                supported_results[prop] = result
        calc = SinglePointDFTCalculator(atoms, **supported_results)

        # This is important! We want to make sure doing something like
        # atoms*(2,2,2) throws away the prior calculator results
        # otherwise we can't do things like run a new calculation
        # with atoms.get_potential_energy() after the transformation
        calc.discard_results_on_any_change = True

        atoms.calc = calc
    return atoms


def make_conventional_cell(atoms):
    """
    Function to make a conventional cell from an Atoms object.

    Args:
        atoms (ase.Atoms): Atoms object

    Returns:
        conventional_atoms (ase.Atoms): Atoms object with a conventional cell
    """

    if isinstance(atoms, Atoms):
        struct = AseAtomsAdaptor.get_structure(atoms)
    else:
        struct = atoms
    conventional_struct = SpacegroupAnalyzer(
        struct
    ).get_conventional_standard_structure()
    conventional_atoms = AseAtomsAdaptor.get_atoms(conventional_struct)

    return conventional_atoms


def invert_slab(slab, return_struct=False):
    """
    Function to invert a Pymatgen slab object, keeping the vacuum
    space in place.

    Args:
        slab (pymatgen.core.surface.Slab): slab to invert
        return_struct (bool): True if a Pymatgen structure (technically, slab) object
        should be returned; False if an ASE atoms object should be returned
            Defaults to False

    Returns:
        inverted_slab (ase.Atoms or pymatgen.core.surface.Slab): inverted slab

    """
    if isinstance(slab, Slab):
        slab_struct = slab
    else:
        raise TypeError("slab must be a pymatgen.core.surface.Slab object")
    frac_coords = slab_struct.frac_coords
    max_z = np.max(frac_coords[:, -1])
    min_z = np.min(frac_coords[:, -1])
    frac_coords[:, -1] = max_z + min_z - frac_coords[:, -1]
    oriented_cell = slab_struct.oriented_unit_cell
    oriented_frac_coords = oriented_cell.frac_coords
    max_oriented_c = np.max(oriented_frac_coords[:, -1])
    min_oriented_c = np.min(oriented_frac_coords[:, -1])
    oriented_frac_coords[:, -1] = (
        max_oriented_c + min_oriented_c - oriented_frac_coords[:, -1]
    )
    inverted_oriented_cell = Structure(
        oriented_cell.lattice,
        oriented_cell.species_and_occu,
        oriented_frac_coords,
    )
    inverted_slab_struct = Slab(
        slab_struct.lattice,
        species=slab_struct.species_and_occu,
        coords=frac_coords,
        miller_index=slab_struct.miller_index,
        oriented_unit_cell=inverted_oriented_cell,
        shift=-slab_struct.shift,
        scale_factor=slab_struct.scale_factor,
    )
    if return_struct:
        inverted_slab = inverted_slab_struct
    else:
        inverted_slab = AseAtomsAdaptor.get_atoms(inverted_slab_struct)

    return inverted_slab


def make_slabs_from_bulk(
    atoms,
    max_index=1,
    min_slab_size=12.0,
    min_vacuum_size=15.0,
    min_length_width=10.0,
    z_fix=2.0,
    return_struct=False,
):
    """
    Function to make slabs from a bulk atoms object.

    Args:
        atoms (ase.Atoms/pymatgen.core.Structre): bulk atoms/structure
        max_index (int): maximum Miller index for slab generation
            Defaults to 1.
        min_slab_size (float): minimum slab size in angstroms
            Defaults to 12.0
        min_vacuum_size (float): minimum vacuum size in angstroms
            Defaults to 15.0
        min_length_width (float): minimum length and width of the slab in angstroms
            Defaults to 10.0
        z_fix (float): distance (in angstroms) from top of slab for which atoms should be fixed
            Defaults to 2.0
        return_struct (bool): True if a Pymatgen structure (technically, slab) object
        should be returned; False if an ASE atoms object should be returned
            Defaults to False

    Returns:
        final_slabs (ase.Atoms or pymatgen.core.surface.Slab): inverted slab
    """

    # Note: This will not work as expected if the slab crosses the
    # unit cell boundary or for 2D systems. See Oxana/Martin's code
    # for the 2D workflow: https://github.com/oxana-a/atomate/blob/ads_wf/atomate/vasp/firetasks/adsorption_tasks.py

    # Use pymatgen to generate slabs
    if isinstance(atoms, Atoms):
        struct = AseAtomsAdaptor.get_structure(atoms)
    else:
        struct = atoms

    slabs = [
        slab
        for slab in generate_all_slabs(
            struct, max_index, min_slab_size, min_vacuum_size
        )
    ]

    # If the two terminations are not equivalent, make new slab
    # by inverting the original slab and add it to the list
    new_slabs = []
    for slab in slabs:
        if not slab.is_symmetric():
            new_slab = invert_slab(slab, return_struct=True)
            new_slabs.append(new_slab)

    slabs.extend(new_slabs)

    # For each slab, make sure the lengths and widths are large enough
    # and fix atoms z_fix away from the top of the slab.
    final_slabs = []
    for slab in slabs:

        # Supercell creation (if necessary)
        a_factor = int(np.ceil(min_length_width / slab.lattice.abc[0]))
        b_factor = int(np.ceil(min_length_width / slab.lattice.abc[1]))
        final_slab = slab * (a_factor, b_factor, 1)

        # Apply constraints by distance from top surface
        # This does not actually create an adsorbate. It is just a
        # useful function for finding surface vs. subsurface sites
        # since you can't just do z_max - z_fix
        if z_fix:
            final_slab = AdsorbateSiteFinder(
                final_slab, selective_dynamics=True, height=z_fix
            ).slab
        final_slabs.append(final_slab)

    if return_struct is False:
        final_slabs = [AseAtomsAdaptor().get_atoms(slab) for slab in final_slabs]

    return final_slabs
