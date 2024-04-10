from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_numbers
from ase.io import read, write
from ase.units import Bohr

from chemsh import *
from chemsh.io.tools import convert_atoms_to_frag

from monty.io import zopen
from monty.os.path import zpath

if TYPE_CHECKING:
    from numpy.typing import NDArray

def get_cluster_info_from_slab(
        adsorbate_slab_file: Path,
        slab_center_idx: list[int],
        adsorbate_idx: list[int] 
) -> tuple[Atoms, Atoms, int, NDArray, NDArray]:
    """
    Read the file containing the periodic slab and adsorbate (geometry optimized) and return the key information needed to create an embedded cluster in ChemShell.

    Parameters
    ----------
    adsorbate_slab_file
        The path to the file containing the adsorbate molecule on the surface slab. It can be in any format that ASE can read.
    adsorbate_idx
        The indices of the atoms that make up the adsorbate molecule.
    slab_center_idx
        The the adsorbate_slab_file, this gives the indices of the atoms that are at the 'center' of the slab right beneath the adsorbate.

    Returns
    -------
    adsorbate
        The Atoms object of the adsorbate molecule.
    slab
        The Atoms object of the surface slab.
    center_position
        The position of the center of the cluster.
    list[float, float, float]
        The vector from the center of the slab to the center of mass of the adsorbate.
    """
    # """
    # Get the necessary information for the cluster from a provided slab file (in any format that ASE can read)
    # """

    adsorbate_slab = read(adsorbate_slab_file)
    # Get the slab 

    # Find indices (within adsorbate_slab) of the slab
    slab_idx = [x for x in list(range(len(adsorbate_slab))) if x not in adsorbate_idx]
    # Create slab from adsorbate_slab
    slab = adsorbate_slab[slab_idx].copy()

    # Find index of the first center atom of the slab as listed in slab_center_idx
    slab_first_atom_idx = [index for index,x in enumerate(slab_idx) if x == slab_center_idx[0]][0]

    # Get the centre of the cluster from the atom indices
    slab_center_position = np.zeros(3)
    for atom_idx in slab_center_idx:
        slab_center_position += adsorbate_slab.get_positions()[atom_idx]

    slab_center_position = slab_center_position/len(slab_center_idx)

    adsorbate = adsorbate_slab[adsorbate_idx].copy()

    # Get the relative distance of the adsorbate from the first center atom of the slab as defined in the slab_center_idx
    adsorbate_com = adsorbate.get_center_of_mass()
    adsorbate_vector_from_slab = adsorbate[0].position - adsorbate_slab[slab_center_idx[0]].position

    # Add the height of the adsorbate from the slab along the z-direction relative to the first center atom of the slab as defined in the slab_center_idx
    adsorbate_com_z_disp = adsorbate_com[2] - adsorbate_slab[slab_center_idx[0]].position[2]
    center_position = np.array([0.0, 0.0, adsorbate_com_z_disp]) +  slab_center_position - adsorbate_slab[slab_center_idx[0]].position

    return adsorbate, slab, slab_first_atom_idx, center_position, adsorbate_vector_from_slab

def run_chemshell(
        slab: Atoms,
        slab_center_idx: int,
        atom_oxi_states: dict[str, float],
        chemsh_radius_active: float=40.0,
        chemsh_radius_cluster: float=60.0,
        chemsh_bq_layer: float=6.0,
        write_xyz_file: bool = False
):
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
    chemsh_radius_cluster
        The radius of the total embedded cluster in Angstroms.
    chemsh_radius_active
        The radius of the active region in Angstroms. This 'active' region is simply region where the charge fitting is performed to ensure correct Madelung potential; it can be a relatively large value.
    chemsh_bq_layer
        The height above the surface to place some additional fitting point charges in Angstroms; simply for better reproduction of the electrostatic potential close to the adsorbate.
    write_xyz_file
        Whether to write an XYZ file of the cluster for visualisation.
    """
    
    # Translate slab such that first Mg atom is at 0,0,0
    slab.translate(-slab.get_positions()[slab_center_idx])

    # Convert ASE Atoms to ChemShell Fragment object
    slab_frag = convert_atoms_to_frag(slab, connect_mode='ionic', dim='2D')

    # Add the atomic charges to the fragment
    slab_frag.addCharges(atom_oxi_states)
    
    # Create the chemshell cluster (i.e., add electrostatic fitting charges) from the fragment
    chemsh_embedded_cluster = slab_frag.construct_cluster(origin=slab_center_idx, 
                                         radius_cluster=chemsh_radius_cluster/Bohr, 
                                         radius_active=chemsh_radius_active/Bohr, 
                                         bq_layer=chemsh_bq_layer/Bohr,
                                         adjust_charge='coordination_scaled',)

    # Save the final cluster to a .pun file
    chemsh_embedded_cluster.save('ChemShell_cluster.pun', 'pun')

    if write_xyz_file:
        # XYZ for visualisation
        chemsh_embedded_cluster.save('ChemShell_cluster.xyz', 'xyz')


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
    Returns a list of list containing the indices of the atoms (in embedded_cluster) which form the quantum clusters in the SKZCAM protocol. The number of clusters created is controlled by the rdf_max parameter.

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
    quantum_cluster_idx = []
    dummy_cation_idx = []
    dummy_anion_idx = []
    for shell_idx in range(shell_max):
        dummy_cation_idx += cation_shells_idx[shell_idx]
        dummy_anion_idx += anion_coord_idx[shell_idx]
        quantum_cluster_idx += [list(set(dummy_cation_idx + dummy_anion_idx))]

    # Get the ECP region for each quantum cluster
    ecp_region_idx = _get_ecp_region(
        embedded_cluster, quantum_cluster_idx, embedded_cluster_all_dist,  ecp_dist
    )

    # Write the quantum clusters to files
    if write_clusters:
        for idx, cluster in enumerate(quantum_cluster_idx):
            cluster_atoms = embedded_cluster[cluster]
            if write_include_ecp:
                ecp_atoms = embedded_cluster[ecp_region_idx[idx]]
                ecp_atoms.set_chemical_symbols(np.array(["U"] * len(ecp_atoms)))
                cluster_atoms = cluster_atoms.copy() + ecp_atoms.copy()
            write(
                Path(write_clusters_path, f"SKZCAM_cluster_{idx}.xyz"),
                cluster_atoms,
            )

    return embedded_cluster, quantum_cluster_idx, ecp_region_idx


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
        raw_pun_file = [line.rstrip().decode("utf-8") if isinstance(line, bytes) else line.rstrip() for line in f]

    # Get the number of atoms and number of atomic charges in the .pun file
    n_atoms = int(raw_pun_file[3].split()[-1])
    n_charges = int(raw_pun_file[4 + n_atoms - 1 + 3].split()[-1])

    # Check if number of atom charges same as number of atom positions
    if n_atoms != n_charges:
        raise ValueError(
            "Number of atomic positions and atomic charges in the .pun file are not the same."
        )

    raw_atom_positions = raw_pun_file[4 : 4 + n_atoms]
    raw_charge_list = raw_pun_file[7 + n_atoms : 7 + 2 * n_atoms]
    charge_list = [float(charge) for charge in raw_charge_list]

    # Add the atomic positions the embedded_cluster Atoms object (converting from Bohr to Angstrom)
    atom_type_list = []
    atom_type_list = []
    atom_number_list = []
    atom_position_list = []
    # Add the atomic positions the embedded_cluster Atoms object (converting from Bohr to Angstrom)
    for _, line in enumerate(raw_atom_positions):
        line_info = line.split()

        # Add the atom type to the atom_type_list
        if line_info[0] in atom_type_dict:
            atom_type_list.append(atom_type_dict[line_info[0]])
        elif line_info[0] == "F":
            atom_type_list.append("fitting charge")
        else:
            atom_type_list.append("unknown")

        # Add the atom number to the atom_number_list and position to the atom_position_list
        atom_number_list += [atomic_numbers[line_info[0]]]
        atom_position_list += [
            [
                float(line_info[1]) * Bohr,
                float(line_info[2]) * Bohr,
                float(line_info[3]) * Bohr,
            ]
        ]

    embedded_cluster = Atoms(numbers=atom_number_list, positions=atom_position_list)

    # Center the embedded cluster so that atom index 0 is at the [0, 0, 0] position
    embedded_cluster.translate(-embedded_cluster[0].position)

    # Add the `oxi_states` and `atom_type` arrays to the Atoms object
    embedded_cluster.set_array("oxi_states", np.array(charge_list))
    embedded_cluster.set_array("atom_type", np.array(atom_type_list))

    return embedded_cluster


def _get_atom_distances(embedded_cluster: Atoms, center_position: NDArray) -> NDArray:
    """
    Returns the distance of all atoms from the center position of the embedded cluster

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates of the embedded cluster.

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
    shells_idx = []

    # Sort the points by distance from the cluster center for the cations only
    distances_sorted = []
    distances_sorted_idx = []
    for i in np.argsort(distances):
        if embedded_cluster.get_array("atom_type")[i] == "cation":
            distances_sorted.append(distances[i])
            distances_sorted_idx.append(i)

    curr_point = distances_sorted[0]
    curr_shell = [curr_point]
    curr_shell_idx = [distances_sorted_idx[0]]

    for idx, point in enumerate(distances_sorted[1:]):
        if point <= curr_point + shell_width:
            curr_shell.append(point)
            curr_shell_idx.append(distances_sorted_idx[idx + 1])
        else:
            shells.append(curr_shell)
            shells_idx.append(curr_shell_idx)
            curr_shell = [point]
            curr_shell_idx = [distances_sorted_idx[idx + 1]]
        curr_point = point
    shells.append(curr_shell)
    shells_idx.append(curr_shell_idx)

    return shells, shells_idx


def _get_anion_coordination(
    embedded_cluster: Atoms,
    cation_shell_idx: list[int],
    dist_matrix: NDArray,
    bond_dist: float = 2.5,
) -> list[int]:
    """
    Returns a list of lists containing the indices of the anions coordinating the cation indices provided.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    cation_shell_idx
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
    anion_coord_idx = []

    # Iterate over the cation shell indices and find the atoms within the bond distance of each cation
    for atom_idx in cation_shell_idx:
        anion_coord_idx += [
            idx
            for idx, dist in enumerate(dist_matrix[atom_idx])
            if (
                dist < bond_dist
                and embedded_cluster.get_array("atom_type")[idx] == "anion"
            )
        ]

    return list(set(anion_coord_idx))


def _get_ecp_region(
    embedded_cluster: Atoms,
    quantum_cluster_idx: list[int],
    dist_matrix: NDArray,
    ecp_dist: float = 6.0,
) -> list[list[int]]:
    """
    Returns a list of lists containing the indices of the atoms in the ECP region of the embedded cluster for each quantum cluster

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    quantum_cluster_idx
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

    ecp_region_idx = []
    dummy_cation_idx = []

    # Iterate over the quantum clusters and find the atoms within the ECP distance of each quantum cluster
    for cluster in quantum_cluster_idx:
        dummy_cation_idx += cluster
        cluster_ecp_region_idx = []
        for atom_idx in dummy_cation_idx:
            for idx, dist in enumerate(dist_matrix[atom_idx]):
                # Check if the atom is within the ecp_dist region and is not in the quantum cluster and is a cation
                if (
                    dist < ecp_dist
                    and idx not in dummy_cation_idx
                    and embedded_cluster.get_array("atom_type")[idx] == "cation"
                ):
                    cluster_ecp_region_idx += [idx]

        ecp_region_idx += [list(set(cluster_ecp_region_idx))]

    return ecp_region_idx
