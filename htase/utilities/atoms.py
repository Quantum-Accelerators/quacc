from ase.atom import Atom
from pymatgen.io.ase import AseAtomsAdaptor
from ase.constraints import FixAtoms
from pymatgen.core.surface import generate_all_slabs, Slab
from pymatgen.core import Structure
from ase.atoms import Atoms
import numpy as np


def invert_slab(slab, return_atoms=True):
    if isinstance(slab, Atoms):
        slab_struct = AseAtomsAdaptor.get_structure(slab)
    else:
        slab_struct = slab
    frac_coords = slab_struct.frac_coords
    max_c = np.max(frac_coords[:, -1])
    min_c = np.min(frac_coords[:, -1])
    frac_coords[:, -1] = max_c + min_c - frac_coords[:, -1]
    oriented_cell = slab_struct.oriented_unit_cell
    oriented_frac_coords = oriented_cell.frac_coords
    max_oriented_c = np.max(oriented_frac_coords[:, -1])
    min_oriented_c = np.min(oriented_frac_coords[:, -1])
    oriented_frac_coords[:, -1] = (
        max_oriented_c + min_oriented_c - oriented_frac_coords[:, -1]
    )
    inverted_oriented_cell = Structure(
        oriented_cell.lattice, oriented_cell.species_and_occu, oriented_frac_coords,
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
    if return_atoms:
        inverted_slab = AseAtomsAdaptor.get_atoms(inverted_slab_struct)
    return inverted_slab


def make_slabs_from_bulk(
    atoms,
    max_index=1,
    min_slab_size=12.0,
    min_vacuum_size=20.0,
    min_length_width=10.0,
    z_fix=2.0,
):
    """
    Function to make slabs from a bulk atoms object.

    Args:
        atoms (ase.Atoms): bulk atoms object
        max_index (int): maximum Miller index for slab generation
            Defaults to 1.
        min_slab_size (float): minimum slab size in angstroms
            Defaults to 12.0
        min_vacuum_size (float): minimum vacuum size in angstroms
            Defaults to 20.0
        min_length_width (float): minimum length and width of the slab in angstroms
            Defaults to 10.0
        z_fix (float): distance (in angstroms) from top of slab for which atoms should be fixed
            Defaults to 2.0

    Returns:
        slabs (list of ase.Atoms): list of slabs
    """

    # Note: This will not work as expected if the slab crosses the
    # unit cell boundary.

    # Use pymatgen to generate slabs
    struct = AseAtomsAdaptor.get_structure(atoms)
    slab_structs = [
        slab_struct
        for slab_struct in generate_all_slabs(
            struct, max_index, min_slab_size, min_vacuum_size
        )
    ]

    # If the two terminations are not equivalent, make new slab
    # by inverting the original slab and add it to the list
    new_slab_structs = []
    for slab_struct in slab_structs:
        if not slab_struct.is_symmetric():
            new_slab_struct = invert_slab(slab_struct, return_atoms=False)
            new_slab_structs.append(new_slab_struct)

    slab_structs.extend(new_slab_structs)
    slabs = [AseAtomsAdaptor().get_atoms(slab_struct) for slab_struct in slab_structs]

    # For each slab, make sure the lengths and widths are large enough
    # and fix atoms z_fix away from the top of the slab.
    for slab in slabs:

        # Supercell creation (if necessary)
        a_factor = np.ceil(min_length_width / slab.cell.lengths[0])
        b_factor = np.ceil(min_length_width / slab.cell.lengths[1])
        slab *= (a_factor, b_factor, 1)

        # Apply constraints by distance from surface
        max_z = np.max(atoms.positions[:, -1])
        constraints = FixAtoms(
            indices=[atom.index for atom in atoms if atom.z < max_z - z_fix]
        )
        slab.set_constraint(constraints)

    return slabs
