import warnings
from copy import deepcopy
from typing import Any, Dict, List

import numpy as np
from ase.atoms import Atoms
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core import Structure
from pymatgen.core.surface import Slab, center_slab, generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor

# NOTES:
# - Anytime an Atoms object is converted to a pmg structure, make sure
# to reattach any .info flags to the Atoms object, e.g. via `new_atoms.info = atoms.info.copy()``.
# Note that atoms.info is mutable, so copy it!
# - All major functions should take in Atoms by default and reutrn Atoms
# by default. Pymatgen structures can be returned with an optional kwarg.
# - If you modify the properties of an input Atoms object in any way, make sure to do so
# on a deepcopy because Atoms objects are mutable.


def flip_atoms(
    atoms: Atoms | Structure | Slab, return_struct: bool = False
) -> Atoms | Structure | Slab:
    """
    Convenience function for vertically flipping periodic atoms or structures

    Parameters
    ----------
    atoms
        .Atoms/.Structure to flip
    return_struct
        True if a Pymatgen structure object should be returned.
        False if an ASE atoms object should be returned

    Returns
    -------
    Atoms | Structure | Slab
        Inverted slab
    """

    if isinstance(atoms, Atoms):
        new_atoms = deepcopy(atoms)
        atoms_info = atoms.info.copy()
    else:
        new_atoms = AseAtomsAdaptor.get_atoms(atoms)
        atoms_info = {}

    new_atoms.rotate(180, "x")
    new_atoms.wrap()

    new_atoms.info = atoms_info
    if return_struct:
        new_atoms = AseAtomsAdaptor.get_structure(new_atoms)

    return new_atoms


def make_slabs_from_bulk(
    atoms: Atoms,
    max_index: int = 1,
    min_slab_size: float = 10.0,
    min_length_width: float = 8.0,
    min_vacuum_size: float = 20.0,
    z_fix: None | float = 2.0,
    flip_asymmetric: bool = True,
    allowed_surface_symbols: List[str] | str = None,
    **slabgen_kwargs,
) -> List[Atoms]:
    """
    Function to make slabs from a bulk atoms object.

    Parameters
    ----------
    atoms
        bulk atoms
    max_index
        Maximum Miller index for slab generation
    min_slab_size
        Minimum slab size (depth) in angstroms
    min_length_width
        Minimum length and width of the slab in angstroms
    min_vacuum_size
        Minimum vacuum size in angstroms
    z_fix
        Distance (in angstroms) from top of slab for which atoms should be fixed
    flip_asymmetric
        If an asymmetric surface should be flipped and added to the list
    allowed_surface_symbols
        List of chemical symbols that must be present on the surface of the slab otherwise the slab will be discarded, e.g. ["Cu", "Ni"]
    **slabgen_kwargs: keyword arguments to pass to the pymatgen generate_all_slabs() function

    Returns
    -------
    List[.Atoms]
        All generated slabs
    """

    # Note: This will not work properly for 2D structures. See Oxana/Martin's code
    # for adjustments for 2D: https://github.com/oxana-a/atomate/blob/ads_wf/atomate/vasp/firetasks/adsorption_tasks.py

    # Use pymatgen to generate slabs
    struct = AseAtomsAdaptor.get_structure(atoms)
    atoms_info = atoms.info.copy()

    if isinstance(allowed_surface_symbols, str):
        allowed_surface_symbols = [allowed_surface_symbols]

    # Make all the slabs
    slabs = generate_all_slabs(
        struct,
        max_index,
        min_slab_size,
        min_vacuum_size,
        center_slab=True,
        **slabgen_kwargs,
    )

    # If the two terminations are not equivalent, make new slab
    # by inverting the original slab and add it to the list
    if flip_asymmetric:

        new_slabs = []
        for slab in slabs:
            if not slab.is_symmetric():

                # Flip the slab and its oriented unit cell
                new_slab = flip_atoms(slab, return_struct=True)
                new_oriented_unit_cell = flip_atoms(
                    slab.oriented_unit_cell, return_struct=True
                )

                # Reconstruct the full slab object, noting the new
                # shift and oriented unit cell
                new_slab = Slab(
                    new_slab.lattice,
                    new_slab.species,
                    coords=new_slab.frac_coords,
                    miller_index=slab.miller_index,
                    oriented_unit_cell=new_oriented_unit_cell,
                    shift=-slab.shift,
                    scale_factor=slab.scale_factor,
                    site_properties=new_slab.site_properties,
                )

                # It looks better to center the inverted slab so we do
                # that here.
                new_slab = center_slab(new_slab)

                # Add the new slab to the list
                new_slabs.append(new_slab)

        slabs.extend(new_slabs)

    # For each slab, make sure the lengths and widths are large enough
    # and fix atoms z_fix away from the top of the slab.
    slabs_with_props = []
    for slab in slabs:

        # Supercell creation (if necessary)
        a_factor = int(np.ceil(min_length_width / slab.lattice.abc[0]))
        b_factor = int(np.ceil(min_length_width / slab.lattice.abc[1]))
        slab_with_props = slab.copy()
        slab_with_props.make_supercell([a_factor, b_factor, 1])

        # Apply constraints by distance from top surface
        # This does not actually create an adsorbate. It is just a
        # useful function for finding surface vs. subsurface sites
        if z_fix:
            slab_with_props = AdsorbateSiteFinder(
                slab_with_props, selective_dynamics=True, height=z_fix
            ).slab

            surface_species = [
                site.specie.symbol
                for site in slab_with_props
                if site.properties["surface_properties"] == "surface"
            ]

            # Check that the desired atoms are on the surface
            if allowed_surface_symbols and not any(
                allowed_surface_atom in surface_species
                for allowed_surface_atom in allowed_surface_symbols
            ):
                continue

        # Add slab to list
        slabs_with_props.append(slab_with_props)

    final_slabs = []
    if len(slabs_with_props) == 0:
        return final_slabs

    # Make atoms objects and store slab stats
    for slab_with_props in slabs_with_props:
        final_slab = AseAtomsAdaptor.get_atoms(slab_with_props)
        slab_stats = {
            "bulk": atoms,
            "miller_index": slab_with_props.miller_index,
            "shift": round(slab_with_props.shift, 3),
            "scale_factor": slab_with_props.scale_factor,
        }
        final_slab.info = atoms_info.copy()
        final_slab.info["slab_stats"] = slab_stats
        final_slabs.append(final_slab)

    return final_slabs


def make_max_slabs_from_bulk(
    atoms: Atoms,
    max_slabs: int = None,
    max_index: int = 1,
    min_slab_size: float = 10.0,
    min_length_width: float = 8.0,
    min_vacuum_size: float = 20.0,
    z_fix: float = 2.0,
    flip_asymmetric: bool = True,
    allowed_surface_symbols: List[str] | str = None,
    **slabgen_kwargs,
) -> List[Atoms]:

    """
    Generate no more than max_slabs number of slabs from a bulk structure.
    The procedure is as follows:
    1. Generate all slabs
    2. If number of slabs is greater than max_slabs, tune ftol from 0.1 to 0.8
    in increments of 0.1. This reduces the number of vertical shifts to consider.
    3. If number of slabs is still greater than max_slabs, only return the slabs
    with the fewest number of atoms per cell such that the returned amount is
    less than or equal to max_slabs.

    Parameters
    ----------
    atoms
        Bulk structure to generate slabs from
    max_slabs
        Maximum number of slabs to generate
    max_index
        Maximum Miller index for slab generation
    min_slab_size
        Minimum slab size (depth) in angstroms
    min_length_width
        Minimum length and width of the slab in angstroms
    min_vacuum_size
        Minimum vacuum size in angstroms
    z_fix
        Distance (in angstroms) from top of slab for which atoms should be fixed
    flip_asymmetric
        If an asymmetric surface should be flipped and added to the list
    allowed_surface_symbols
        List of chemical symbols that must be present on the surface of the slab otherwise
        the slab will be discarded, e.g. ["Cu", "Ni"]
    **slabgen_kwargs: keyword arguments to pass to the pymatgen generate_all_slabs() function

    Returns:
    --------
    List[.Atoms]
        List of slabs

    """

    if isinstance(allowed_surface_symbols, str):
        allowed_surface_symbols = [allowed_surface_symbols]

    slabs = make_slabs_from_bulk(
        atoms,
        max_index=max_index,
        min_slab_size=min_slab_size,
        min_length_width=min_length_width,
        min_vacuum_size=min_vacuum_size,
        z_fix=z_fix,
        flip_asymmetric=flip_asymmetric,
        allowed_surface_symbols=allowed_surface_symbols,
        **slabgen_kwargs,
    )

    # Try to reduce the number of slabs if the user really wants it...
    # (desperate times call for desperate measures)
    if max_slabs and slabs is not None and len(slabs) > max_slabs:

        if len(slabs) > max_slabs:
            warnings.warn(
                f"You requested {max_slabs} slabs, but {len(slabs)} were generated. Tuning ftol in generate_all_slabs() to try to reduce the number of slabs, at the expense of sampling fewer surface configurations.",
            )
            for ftol in np.arange(0.1, 0.9, 0.1):
                slabgen_kwargs["ftol"] = ftol
                slabs_ftol = make_slabs_from_bulk(
                    atoms,
                    max_index=max_index,
                    min_slab_size=min_slab_size,
                    min_length_width=min_length_width,
                    min_vacuum_size=min_vacuum_size,
                    z_fix=z_fix,
                    flip_asymmetric=flip_asymmetric,
                    allowed_surface_symbols=allowed_surface_symbols,
                    **slabgen_kwargs,
                )
                if len(slabs_ftol) < len(slabs):
                    slabs = slabs_ftol
                if len(slabs) <= max_slabs:
                    break

        if len(slabs) > max_slabs:
            warnings.warn(
                f"You requested {max_slabs} slabs, but {len(slabs)} were generated. Could not reduce further. Picking the smallest slabs by number of atoms.",
            )
            slabs.sort(key=len)
            slabs = slabs[0:max_slabs]

    return slabs


def make_adsorbate_structures(
    atoms: Atoms,
    adsorbate: Atoms,
    min_distance: float = 2.0,
    modes: List[str] | str = None,
    allowed_surface_symbols: List[str] | str = None,
    allowed_surface_indices: List[int] | int = None,
    ads_site_finder_kwargs: Dict[str, Any] = None,
    find_ads_sites_kwargs: Dict[str, Any] = None,
) -> List[Atoms]:
    """
    Add a single adsorbate to a structure for every requested adsorption mode

    Parameters
    ----------
    atoms
        The atoms to add adsorbates to.
    adsorbate
        The adsorbate to add.
        Note: It will be placed on the surface in the exact input orientation provided by the user (the adsorption mode is
        along the c axis and the coordinating atom is the one in the -z direction).
    min_distance
        The distance between the adsorbate and the surface site.
    modes
        The adsorption mode(s) to consider. Options include: "ontop", "bridge", "hollow", "subsurface".
    allowed_surface_symbols
        The symbols of surface atoms to consider. If None, will use all surface atoms.
    allowed_surface_indices
        The indices of surface atoms to consider. If None, will use all surface atoms. Generally used if a specific site is to be excluded from the set.
    ads_site_finder_kwargs
        The keyword arguments to pass to the AdsorbateSiteFinder().
    find_ads_sites_kwargs
        The keyword arguments to pass to AdsorbateSiteFinder.find_adsorption_sites().

    Returns:
    --------
    List[ase.Atoms]
        The structures with adsorbates

    """
    if modes is None:
        modes = ["ontop", "bridge", "hollow"]
    if isinstance(modes, str):
        modes = [modes]
    if isinstance(allowed_surface_symbols, str):
        allowed_surface_symbols = [allowed_surface_symbols]
    if isinstance(allowed_surface_indices, int):
        allowed_surface_indices = [allowed_surface_indices]

    ads_site_finder_kwargs = ads_site_finder_kwargs or {}
    find_ads_sites_kwargs = find_ads_sites_kwargs or {}

    # Check for double-used parameters
    if min_distance and "distance" in find_ads_sites_kwargs:
        raise ValueError(
            "Cannot specify both min_distance and find_ads_sites_kwargs['distance']",
        )
    if modes and "positions" in find_ads_sites_kwargs:
        raise ValueError(
            "Cannot specify both modes and find_ads_sites_kwargs['positions']",
        )
    find_ads_sites_kwargs["distance"] = min_distance
    find_ads_sites_kwargs["positions"] = [mode.lower() for mode in modes]

    # Check the provided surface indices are reasonable
    atom_indices = [atom.index for atom in atoms]
    if allowed_surface_indices and not all(
        idx in atom_indices for idx in allowed_surface_indices
    ):
        raise ValueError(
            "All indices in allowed_surface_indices must be in atoms.",
            allowed_surface_indices,
            atom_indices,
        )

    # Add 0.0 initial magmoms to atoms/adsorbate if needed
    if atoms.has("initial_magmoms") and not adsorbate.has("initial_magmoms"):
        adsorbate.set_initial_magnetic_moments([0.0] * len(adsorbate))
    if adsorbate.has("initial_magmoms") and not atoms.has("initial_magmoms"):
        atoms.set_initial_magnetic_moments([0.0] * len(atoms))

    # Make a Pymatgen structure and molecule
    struct = AseAtomsAdaptor.get_structure(atoms)
    mol = AseAtomsAdaptor.get_molecule(adsorbate, charge_spin_check=False)

    # Get the adsorption sites
    ads_finder = AdsorbateSiteFinder(struct, **ads_site_finder_kwargs)
    ads_sites = ads_finder.find_adsorption_sites(**find_ads_sites_kwargs)

    # Find and add the adsorbates
    new_atoms = []
    for mode, ads_coords in ads_sites.items():

        # Check if mode is in desired list
        if mode not in modes:
            continue

        for ads_coord in ads_coords:

            # Place adsorbate
            struct_with_adsorbate = ads_finder.add_adsorbate(mol, ads_coord)

            # Convert back to Atoms object
            atoms_with_adsorbate = AseAtomsAdaptor.get_atoms(struct_with_adsorbate)

            # Get distance matrix between adsorbate binding atom and surface
            adsorbate_index = len(atoms) + np.argmin(atom.z for atom in adsorbate)
            d = atoms_with_adsorbate.get_all_distances(mic=True)
            d = d[atom_indices, :]
            d = d[:, adsorbate_index]

            # Find closest surface atoms
            min_d = min(d)
            surface_atom_indices = [
                i for i, val in enumerate(d) if (min_d - 0.01) <= val <= (min_d + 0.01)
            ]
            surface_atom_symbols = atoms_with_adsorbate[
                surface_atom_indices
            ].get_chemical_symbols()

            # Check if surface binding site is not in the specified
            # user list. If so, skip this one
            if allowed_surface_symbols and not any(
                surface_atom_symbol in allowed_surface_symbols
                for surface_atom_symbol in surface_atom_symbols
            ):
                continue

            if allowed_surface_indices and not any(
                surface_atom_idx in allowed_surface_indices
                for surface_atom_idx in surface_atom_indices
            ):
                continue

            # Store adsorbate info
            atoms_with_adsorbate.info = atoms.info.copy()
            ads_stats = {
                "adsorbate": adsorbate,
                "initial_mode": mode,
                "surface_atoms_symbols": surface_atom_symbols,
                "surface_atoms_indices": surface_atom_indices,
            }
            if atoms_with_adsorbate.info.get("adsorbates", None) is None:
                atoms_with_adsorbate.info["adsorbates"] = [ads_stats]
            else:
                atoms_with_adsorbate.info["adsorbates"].extend([ads_stats])

            # Add slab+adsorbate to list
            new_atoms.append(atoms_with_adsorbate)

    return new_atoms


def get_surface_energy(
    bulk: Atoms, slab: Atoms, bulk_energy: float, slab_energy: float
):
    """
    Calculate the surface energy to form a given surface slab from a bulk structure.
    For asymmetric slabs, this is best thought of as an average surface energy, where
    the cleavage energy would be 2x this value.

    Parameters:
    -----------
    bulk
        The bulk structure.
    slab
        The slab structure.
    bulk_energy
        The total energy of the bulk structure.
    slab_energy
        The total energy of the slab structure.
    """

    alpha = len(slab) / len(bulk)
    matrix = slab.cell.array
    A = np.linalg.norm(np.cross(matrix[0, :], matrix[1, :]))
    surface_energy = (slab_energy - alpha * bulk_energy) / (2 * A)

    return surface_energy
