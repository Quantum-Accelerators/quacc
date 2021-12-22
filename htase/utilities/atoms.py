from pymatgen.io.ase import AseAtomsAdaptor
from ase.constraints import FixAtoms
from pymatgen.core.surface import generate_all_slabs
import numpy as np


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

    # Use pymatgen to generate slabs
    struct = AseAtomsAdaptor.get_structure(atoms)
    slabs = [
        AseAtomsAdaptor.get_atoms(slab)
        for slab in generate_all_slabs(
            struct, max_index, min_slab_size, min_vacuum_size
        )
    ]

    # For each slab, make sure the lengths and widths are large enough
    # and fix atoms z_fix away from the top of the slab.
    for slab in slabs:

        # Supercell creation (if necessary)
        a_factor = np.ceil(min_length_width / slab.cell.lengths[0])
        b_factor = np.ceil(min_length_width / slab.cell.lengths[1])
        slab *= (a_factor, b_factor, 1)

        # Apply constraints by distance from surface
        # TODO: This will not work as expected if the slab crosses the
        # unit cell boundary.
        max_z = np.max(atoms.positions[:, -1])
        constraints = FixAtoms(
            indices=[atom.index for atom in atoms if atom.z < max_z - z_fix]
        )
        slab.set_constraint(constraints)

    return slabs
