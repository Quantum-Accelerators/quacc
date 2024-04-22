from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_numbers
from ase.io import read, write
from ase.units import Bohr
from monty.dev import requires
from monty.io import zopen
from monty.os.path import zpath

if TYPE_CHECKING:
    from numpy.typing import NDArray

has_chemshell = find_spec("chemsh") is not None


def get_cluster_info_from_slab(
    adsorbate_slab_file: str | Path,
    slab_center_indices: list[int],
    adsorbate_indices: list[int],
) -> tuple[Atoms, Atoms, int, NDArray, NDArray]:
    """
    Read the file containing the periodic slab and adsorbate (geometry optimized) and return the key information needed to create an embedded cluster in ChemShell.

    Parameters
    ----------
    adsorbate_slab_file
        The path to the file containing the adsorbate molecule on the surface slab. It can be in any format that ASE can read.
    adsorbate_indices
        The indices of the atoms that make up the adsorbate molecule.
    slab_center_indices
        The the adsorbate_slab_file, this gives the indices of the atoms that are at the 'center' of the slab right beneath the adsorbate.

    Returns
    -------
    Atoms
        The Atoms object of the adsorbate molecule.
    Atoms
        The Atoms object of the surface slab.
    int
        The index of the first atom of the slab as listed in slab_center_indices.
    NDArray
        The position of the center of the cluster.
    NDArray
        The vector from the center of the slab to the center of mass of the adsorbate.
    """

    # Get the necessary information for the cluster from a provided slab file (in any format that ASE can read)
    adsorbate_slab = read(adsorbate_slab_file)

    # Find indices (within adsorbate_slab) of the slab
    slab_indices = [
        i for i, _ in enumerate(adsorbate_slab) if i not in adsorbate_indices
    ]

    # Create slab from adsorbate_slab
    slab = adsorbate_slab[slab_indices]

    # Find index of the first center atom of the slab as listed in slab_center_indices
    slab_first_atom_idx = next(
        index for index, x in enumerate(slab_indices) if x == slab_center_indices[0]
    )

    # Get the center of the cluster from the atom indices
    slab_center_position = adsorbate_slab[slab_center_indices].get_positions().sum(
        axis=0
    ) / len(slab_center_indices)

    adsorbate = adsorbate_slab[adsorbate_indices]

    # Get the relative distance of the adsorbate from the first center atom of the slab as defined in the slab_center_indices
    adsorbate_com = adsorbate.get_center_of_mass()
    adsorbate_vector_from_slab = (
        adsorbate[0].position - adsorbate_slab[slab_center_indices[0]].position
    )

    # Add the height of the adsorbate from the slab along the z-direction relative to the first center atom of the slab as defined in the slab_center_indices
    adsorbate_com_z_disp = (
        adsorbate_com[2] - adsorbate_slab[slab_center_indices[0]].position[2]
    )
    center_position = (
        np.array([0.0, 0.0, adsorbate_com_z_disp])
        + slab_center_position
        - adsorbate_slab[slab_center_indices[0]].position
    )

    return (
        adsorbate,
        slab,
        slab_first_atom_idx,
        center_position,
        adsorbate_vector_from_slab,
    )


@requires(has_chemshell, "ChemShell is not installed")
def generate_chemshell_cluster(
    slab: Atoms,
    slab_center_idx: int,
    atom_oxi_states: dict[str, float],
    filepath: str | Path,
    chemsh_radius_active: float = 40.0,
    chemsh_radius_cluster: float = 60.0,
    chemsh_bq_layer: float = 6.0,
    write_xyz_file: bool = False,
) -> None:
    """
    Run ChemShell to create an embedded cluster from a slab.

    Parameters
    ----------
    slab
        The Atoms object of the slab.
    slab_center_idx
        The index of the (first) atom at the center of the slab, this index corresponds to the atom in the slab_center_idx list but adjusted for the slab (which does not contain the adsorbate atoms)
    atom_oxi_states
        The oxidation states of the atoms in the slab as a dictionary
    filepath
        The location where the ChemShell output files will be written.
    chemsh_radius_cluster
        The radius of the total embedded cluster in Angstroms.
    chemsh_radius_active
        The radius of the active region in Angstroms. This 'active' region is simply region where the charge fitting is performed to ensure correct Madelung potential; it can be a relatively large value.
    chemsh_bq_layer
        The height above the surface to place some additional fitting point charges in Angstroms; simply for better reproduction of the electrostatic potential close to the adsorbate.
    write_xyz_file
        Whether to write an XYZ file of the cluster for visualisation.

    Returns
    -------
    None
    """
    from chemsh.io.tools import convert_atoms_to_frag

    # Translate slab such that first Mg atom is at 0,0,0
    slab.translate(-slab.get_positions()[slab_center_idx])

    # Convert ASE Atoms to ChemShell Fragment object
    slab_frag = convert_atoms_to_frag(slab, connect_mode="ionic", dim="2D")

    # Add the atomic charges to the fragment
    slab_frag.addCharges(atom_oxi_states)

    # Create the chemshell cluster (i.e., add electrostatic fitting charges) from the fragment
    chemsh_embedded_cluster = slab_frag.construct_cluster(
        origin=slab_center_idx,
        radius_cluster=chemsh_radius_cluster / Bohr,
        radius_active=chemsh_radius_active / Bohr,
        bq_layer=chemsh_bq_layer / Bohr,
        adjust_charge="coordination_scaled",
    )

    # Save the final cluster to a .pun file
    chemsh_embedded_cluster.save(Path(filepath).with_suffix(".pun"), "pun")

    if write_xyz_file:
        # XYZ for visualisation
        chemsh_embedded_cluster.save(Path(filepath).with_suffix(".xyz"), "xyz")


def create_skzcam_clusters(
    pun_file: str | Path,
    center_position: NDArray,
    atom_oxi_states: dict[str, float],
    shell_max: int = 10,
    shell_width: float = 0.1,
    bond_dist: float = 2.5,
    ecp_dist: float = 6.0,
    write_clusters: bool = False,
    write_clusters_path: str | Path = ".",
    write_include_ecp: bool = False,
) -> tuple[Atoms, list[list[int]], list[list[int]]]:
    """
    From a provided .pun file (generated by ChemShell), this function creates quantum clusters using the SKZCAM protocol. It will return the embedded cluster Atoms object and the indices of the atoms in the quantum clusters and the ECP region. The number of clusters created is controlled by the rdf_max parameter.

    Parameters
    ----------
    pun_file
        The path to the .pun file created by ChemShell to be read.
    center_position
        The position of the center of the embedded cluster (i.e., position of the adsorbate).
    atom_oxi_states
        A dictionary containing the atomic symbols as keys and the oxidation states as values.
    shell_max
        The maximum number of quantum clusters to be created.
    shell_width
        Defines the distance between atoms within shells; this is the maximum distance between any two atoms within the shell.
    bond_dist
        The distance within which an anion is considered to be coordinating a cation.
    ecp_dist
        The distance from edges of the quantum cluster to define the ECP region.
    write_clusters
        If True, the quantum clusters will be written to a file.
    write_clusters_path
        The path to the file where the quantum clusters will be written.
    write_include_ecp
        If True, the ECP region will be included in the quantum clusters.

    Returns
    -------
    Atoms
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file.
    list[list[int]]
        A list of lists containing the indices of the atoms in each quantum cluster.
    list[list[int]]
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster.
    """

    # Read the .pun file and create the embedded_cluster Atoms object
    embedded_cluster = convert_pun_to_atoms(pun_file, atom_oxi_states)

    # Get distances of all atoms from the cluster center
    atom_center_distances = _get_atom_distances(embedded_cluster, center_position)

    # Determine the cation shells from the center of the embedded cluster
    _, cation_shells_idx = _find_cation_shells(
        embedded_cluster, atom_center_distances, shell_width
    )

    # Create the distance matrix for the embedded cluster
    embedded_cluster_all_dist = embedded_cluster.get_all_distances()

    # Create the anion coordination list for each cation shell
    anion_coord_idx = []
    for shell_idx in range(shell_max):
        cation_shell = cation_shells_idx[shell_idx]
        anion_coord_idx += [
            _get_anion_coordination(
                embedded_cluster, cation_shell, embedded_cluster_all_dist, bond_dist
            )
        ]

    # Create the quantum clusters by summing up the indices of the cations and their coordinating anions
    quantum_cluster_indices = []
    dummy_cation_indices = []
    dummy_anion_indices = []
    for shell_idx in range(shell_max):
        dummy_cation_indices += cation_shells_idx[shell_idx]
        dummy_anion_indices += anion_coord_idx[shell_idx]
        quantum_cluster_indices += [
            list(set(dummy_cation_indices + dummy_anion_indices))
        ]

    # Get the ECP region for each quantum cluster
    ecp_region_indices = _get_ecp_region(
        embedded_cluster, quantum_cluster_indices, embedded_cluster_all_dist, ecp_dist
    )

    # Write the quantum clusters to files
    if write_clusters:
        for idx, cluster in enumerate(quantum_cluster_indices):
            cluster_atoms = embedded_cluster[cluster]
            if write_include_ecp:
                ecp_atoms = embedded_cluster[ecp_region_indices[idx]]
                ecp_atoms.set_chemical_symbols(np.array(["U"] * len(ecp_atoms)))
                cluster_atoms = cluster_atoms.copy() + ecp_atoms.copy()
            write(
                Path(write_clusters_path, f"SKZCAM_cluster_{idx}.xyz"),
                cluster_atoms,
            )

    return embedded_cluster, quantum_cluster_indices, ecp_region_indices


def convert_pun_to_atoms(
    pun_file: str | Path, atom_oxi_states: dict[str, float]
) -> Atoms:
    """
    Reads a .pun file and returns an ASE Atoms object containing the atomic coordinates,
    point charges/oxidation states, and atom types.

    Parameters
    ----------
    pun_file
        The path to the .pun file created by ChemShell to be read.
    atom_oxi_states
        A dictionary containing the atomic symbols as keys and the oxidation states as values.

    Returns
    -------
    Atoms
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file.
        The `oxi_states` array contains the atomic charges, and the `atom_type` array contains the
        atom types (cation, anion, neutral).
    """

    # Create an empty embedded_cluster Atoms object
    embedded_cluster = Atoms()

    # Create a dictionary containing the atom types and whether they are cations or anions
    atom_type_dict = {
        atom: "cation" if oxi_state > 0 else "anion" if oxi_state < 0 else "neutral"
        for atom, oxi_state in atom_oxi_states.items()
    }

    # Load the pun file as a list of strings
    with zopen(zpath(Path(pun_file))) as f:
        raw_pun_file = [
            line.rstrip().decode("utf-8") if isinstance(line, bytes) else line.rstrip()
            for line in f
        ]

    # Get the number of atoms and number of atomic charges in the .pun file
    n_atoms = int(raw_pun_file[3].split()[-1])
    n_charges = int(raw_pun_file[4 + n_atoms - 1 + 3].split()[-1])

    # Check if number of atom charges same as number of atom positions
    if n_atoms != n_charges:
        raise ValueError(
            "Number of atomic positions and atomic charges in the .pun file are not the same."
        )

    raw_atom_positions = raw_pun_file[4 : 4 + n_atoms]
    raw_charges = raw_pun_file[7 + n_atoms : 7 + 2 * n_atoms]
    charges = [float(charge) for charge in raw_charges]

    # Add the atomic positions the embedded_cluster Atoms object (converting from Bohr to Angstrom)
    atom_types = []
    atom_numbers = []
    atom_positions = []
    # Add the atomic positions the embedded_cluster Atoms object (converting from Bohr to Angstrom)
    for _, line in enumerate(raw_atom_positions):
        line_info = line.split()

        # Add the atom type to the atom_type_list
        if line_info[0] in atom_type_dict:
            atom_types.append(atom_type_dict[line_info[0]])
        elif line_info[0] == "F":
            atom_types.append("fitting charge")
        else:
            atom_types.append("unknown")

        # Add the atom number to the atom_number_list and position to the atom_position_list
        atom_numbers += [atomic_numbers[line_info[0]]]
        atom_positions += [
            [
                float(line_info[1]) * Bohr,
                float(line_info[2]) * Bohr,
                float(line_info[3]) * Bohr,
            ]
        ]

    embedded_cluster = Atoms(numbers=atom_numbers, positions=atom_positions)

    # Center the embedded cluster so that atom index 0 is at the [0, 0, 0] position
    embedded_cluster.translate(-embedded_cluster[0].position)

    # Add the `oxi_states` and `atom_type` arrays to the Atoms object
    embedded_cluster.set_array("oxi_states", np.array(charges))
    embedded_cluster.set_array("atom_type", np.array(atom_types))

    return embedded_cluster


def insert_adsorbate_to_embedded_cluster(
    embedded_cluster: Atoms,
    adsorbate: Atoms,
    adsorbate_vector_from_slab: NDArray,
    quantum_cluster_indices: list[list[int]] | None = None,
    ecp_region_idx: list[list[int]] | None = None,
) -> tuple[Atoms, list[list[int]], list[list[int]]]:
    """
    Insert the adsorbate into the embedded cluster and update the quantum cluster and ECP region indices.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file.
    adsorbate
        The ASE Atoms object of the adsorbate molecule.
    adsorbate_vector_from_slab
        The vector from the first atom of the embedded cluster to the center of mass of the adsorbate.
    quantum_cluster_indices
        A list of lists containing the indices of the atoms in each quantum cluster.
    ecp_region_idx
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster.

    Returns
    -------
    Atoms
        The ASE Atoms object containing the adsorbate and embedded cluster
    list[list[int]]
        A list of lists containing the indices of the atoms in each quantum cluster.
    list[list[int]]
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster.
    """

    # Remove PBC from the adsorbate
    adsorbate.set_pbc(False)

    # Translate the adsorbate to the correct position relative to the slab
    adsorbate.translate(-adsorbate[0].position + adsorbate_vector_from_slab)

    # Set oxi_state and atom_type arrays for the adsorbate
    adsorbate.set_array("oxi_states", np.array([0.0] * len(adsorbate)))
    adsorbate.set_array("atom_type", np.array(["adsorbate"] * len(adsorbate)))

    # Add the adsorbate to the embedded cluster
    embedded_adsorbate_cluster = adsorbate + embedded_cluster

    # Update the quantum cluster and ECP region indices
    if quantum_cluster_indices is not None:
        quantum_cluster_indices = [
            [idx + len(adsorbate) for idx in cluster]
            for cluster in quantum_cluster_indices
        ]
    if ecp_region_idx is not None:
        ecp_region_idx = [
            [idx + len(adsorbate) for idx in cluster] for cluster in ecp_region_idx
        ]

    return embedded_adsorbate_cluster, quantum_cluster_indices, ecp_region_idx


def _get_atom_distances(embedded_cluster: Atoms, center_position: NDArray) -> NDArray:
    """
    Returns the distance of all atoms from the center position of the embedded cluster

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates of the embedded cluster.
    center_position
        The position of the center of the embedded cluster (i.e., position of the adsorbate).

    Returns
    -------
    NDArray
        An array containing the distances of each atom in the Atoms object from the cluster center.
    """

    return np.array(
        [np.linalg.norm(atom.position - center_position) for atom in embedded_cluster]
    )


def _find_cation_shells(
    embedded_cluster: Atoms,
    distances: NDArray,
    shell_width: float = 0.1,
) -> list[list[int]]:
    """
    Returns a list of lists containing the indices of the cations in each shell, based on distance from the embedded cluster center.
    This is achieved by clustering the data based on the DBSCAN clustering algorithm.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    distances
        The distance of atoms from the cluster center.
    shell_width
        Defines the distance between atoms within shells; this is the maximum distance between any two atoms within the shell

    Returns
    -------
    list[list[int]]
        A list of lists containing the indices of the cations in each shell.
    """

    # Define the empty list to store the cation shells
    shells = []
    shells_indices = []

    # Sort the points by distance from the cluster center for the cations only
    distances_sorted = []
    distances_sorted_indices = []
    for i in np.argsort(distances):
        if embedded_cluster.get_array("atom_type")[i] == "cation":
            distances_sorted.append(distances[i])
            distances_sorted_indices.append(i)

    curr_point = distances_sorted[0]
    curr_shell = [curr_point]
    curr_shell_idx = [distances_sorted_indices[0]]

    for idx, point in enumerate(distances_sorted[1:]):
        if point <= curr_point + shell_width:
            curr_shell.append(point)
            curr_shell_idx.append(distances_sorted_indices[idx + 1])
        else:
            shells.append(curr_shell)
            shells_indices.append(curr_shell_idx)
            curr_shell = [point]
            curr_shell_idx = [distances_sorted_indices[idx + 1]]
        curr_point = point
    shells.append(curr_shell)
    shells_indices.append(curr_shell_idx)

    return shells, shells_indices


def _get_anion_coordination(
    embedded_cluster: Atoms,
    cation_shell_indices: list[int],
    dist_matrix: NDArray,
    bond_dist: float = 2.5,
) -> list[int]:
    """
    Returns a list of lists containing the indices of the anions coordinating the cation indices provided.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    cation_shell_indices
        A list of the indices of the cations in the cluster.
    dist_matrix
        A matrix containing the distances between each pair of atoms in the embedded cluster.
    bond_dist
        The distance within which an anion is considered to be coordinating a cation.

    Returns
    -------
    list[int]
        A list containing the indices of the anions coordinating the cation indices.
    """

    # Define the empty list to store the anion coordination
    anion_coord_indices = []

    # Iterate over the cation shell indices and find the atoms within the bond distance of each cation
    for atom_idx in cation_shell_indices:
        anion_coord_indices += [
            idx
            for idx, dist in enumerate(dist_matrix[atom_idx])
            if (
                dist < bond_dist
                and embedded_cluster.get_array("atom_type")[idx] == "anion"
            )
        ]

    return list(set(anion_coord_indices))


def _get_ecp_region(
    embedded_cluster: Atoms,
    quantum_cluster_indices: list[int],
    dist_matrix: NDArray,
    ecp_dist: float = 6.0,
) -> list[list[int]]:
    """
    Returns a list of lists containing the indices of the atoms in the ECP region of the embedded cluster for each quantum cluster

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    quantum_cluster_indices
        A list of lists containing the indices of the atoms in each quantum cluster.
    dist_matrix
        A matrix containing the distances between each pair of atoms in the embedded cluster.
    ecp_dist
        The distance from edges of the quantum cluster to define the ECP region.

    Returns
    -------
    list[list[int]]
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster.
    """

    ecp_region_indices = []
    dummy_cation_indices = []

    # Iterate over the quantum clusters and find the atoms within the ECP distance of each quantum cluster
    for cluster in quantum_cluster_indices:
        dummy_cation_indices += cluster
        cluster_ecp_region_idx = []
        for atom_idx in dummy_cation_indices:
            for idx, dist in enumerate(dist_matrix[atom_idx]):
                # Check if the atom is within the ecp_dist region and is not in the quantum cluster and is a cation
                if (
                    dist < ecp_dist
                    and idx not in dummy_cation_indices
                    and embedded_cluster.get_array("atom_type")[idx] == "cation"
                ):
                    cluster_ecp_region_idx += [idx]

        ecp_region_indices += [list(set(cluster_ecp_region_idx))]

    return ecp_region_indices
