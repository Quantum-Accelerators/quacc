from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import Atoms
from ase.io import write
from ase.units import Bohr

if TYPE_CHECKING:
    from numpy.typing import NDArray

def convert_pun_to_atoms(
    pun_file: str | Path, atom_oxi_states: dict[str, float]
) -> Atoms:
    """
    Reads a .pun file and returns an ASE Atoms object containing the atomic coordinates, point charges/oxidation states and atom types.

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
    raw_pun_file = [line.rstrip() for line in Path.open(pun_file)]

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
    for _, line in enumerate(raw_atom_positions):
        line_info = line.split()

        # Add the atom type to the atom_type_list
        if line_info[0] in atom_type_dict:
            atom_type_list.append(atom_type_dict[line_info[0]])
        elif line_info[0] == "F":
            atom_type_list.append("fitting charge")
        else:
            atom_type_list.append("unknown")

        atom = Atoms(
            line_info[0],
            positions=[
                [
                    float(line_info[1]) * Bohr,
                    float(line_info[2]) * Bohr,
                    float(line_info[3]) * Bohr,
                ]
            ],
        )
        embedded_cluster = embedded_cluster.copy() + atom.copy()

    # Centre the embedded cluster so that atom index 0 is at the [0, 0, 0] position
    embedded_cluster.translate(-embedded_cluster[0].position)

    # Add the `oxi_states` and `atom_type` arrays to the Atoms object
    embedded_cluster.set_array("oxi_states", np.array(charge_list))
    embedded_cluster.set_array("atom_type", np.array(atom_type_list))

    return embedded_cluster


def get_atom_distances(embedded_cluster: Atoms, center_position: NDArray) -> NDArray:
    """
    Returns the distance of all atoms from the centre position of the embedded cluster

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates of the embedded cluster.

    Returns
    -------
    NDArray
        An array containing the distances of each atom in the Atoms object from the cluster centre.
    """

    return np.array([np.linalg.norm(atom.position - center_position) for atom in embedded_cluster])


def find_cation_shells(
    embedded_cluster: Atoms,
    distances: NDArray,
    shell_width: float = 0.005,
) -> list[list[int]]:
    """
    Returns a list of lists containing the indices of the cations in each shell, based on distance from the embedded cluster centre. This is achieved by clustering the data based on the DBSCAN clustering algorithm.

    Parameters
    ----------
    embedded_cluster : Atoms
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    distances : np.ndarray
        The distance of atoms from the cluster centre.
    shell_width : float
        Defines the distance between atoms within shells; this is the maximum distance between any two atoms within the shell

    Returns
    -------
    shells: list
        A list of lists containing the indices of the cations in each shell.
    """

    # Define the empty list to store the cation shells
    shells = []
    shells_idx = []

    # Sort the points by distance from the cluster centre for the cations only
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


def get_anion_coordination(
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

    Returns
    -------
    anion_coord_idx
        A list containing the indices of the anions coordinating the cation indices.
    """

    # Define the empty list to store the anion coordination
    anion_coord_idx = []

    # Iterate over the cation shell indices and find the atoms within the bond distance of each cation
    for atom_idx in cation_shell_idx:
        anion_coord_idx += [
            idx
            for idx, dist in enumerate(dist_matrix[atom_idx])
            if (dist < bond_dist and embedded_cluster.get_array("atom_type")[idx] == "anion")
        ]

    return list(set(anion_coord_idx))


def get_ecp_region(
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
    ecp_region_idx
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


def create_skzcam_clusters(
    pun_file: str | Path,
    center_position: NDArray,
    atom_oxi_states: dict[str, float],
    shell_max: int = 10,
    shell_width: float = 0.005,
    bond_dist: float = 2.5,
    write_clusters: bool = False,
    write_clusters_path: str | Path | None = None,
    write_cluster_name: str = "SKZCAM_cluster",
    write_include_ecp: bool = False,
) -> list[list[int]]:
    """
    Returns a list of list containing the indices of the atoms (in embedded_cluster) which form the quantum clusters in the SKZCAM protocol. The number of clusters created is controlled by the rdf_max parameter.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    shell_max
        The maximum number of quantum clusters to be created.
    write_clusters
        If True, the quantum clusters will be written to a file.
    write_clusters_path
        The path to the file where the quantum clusters will be written.
        Defaults to "SKZ
    write_cluster_name
        The name of the file where the quantum clusters will be written.
    write_include_ecp
        If True, the ECP region will be included in the quantum clusters.

    Returns
    -------
    clusters
        A list of lists containing the indices of the atoms in each quantum cluster.
    """

    # Read the .pun file and create the embedded_cluster Atoms object
    embedded_cluster = convert_pun_to_atoms(pun_file, atom_oxi_states)

    # Get distances of all atoms from the cluster centre
    atom_centre_distances = get_atom_distances(embedded_cluster, center_position)

    # Determine the cation shells from the centre of the embedded cluster
    _, cation_shells_idx = find_cation_shells(
        embedded_cluster, atom_centre_distances, shell_width
    )

    # Create the distance matrix for the embedded cluster
    embedded_cluster_all_dist = embedded_cluster.get_all_distances()

    anion_coord_idx = []
    # Create the anion coordination list for each cation shell
    for shell_idx in range(shell_max):
        cation_shell = cation_shells_idx[shell_idx]
        anion_coord_idx += [
            get_anion_coordination(embedded_cluster, cation_shell, embedded_cluster_all_dist, bond_dist)
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
    ecp_region_idx = get_ecp_region(embedded_cluster, quantum_cluster_idx, embedded_cluster_all_dist)

    if write_clusters:
        if write_clusters_path is None:
            write_clusters_path = Path.cwd()
        # Write the quantum clusters to a file
        for idx, cluster in enumerate(quantum_cluster_idx):
            cluster_atoms = embedded_cluster[cluster]
            if write_include_ecp:
                ecp_atoms = embedded_cluster[ecp_region_idx[idx]]
                ecp_atoms.set_chemical_symbols(np.array(["U"] * len(ecp_atoms)))
                cluster_atoms = cluster_atoms.copy() + ecp_atoms.copy()

            write(
                Path(f"{write_clusters_path}/{write_cluster_name}_{idx}.xyz"),
                cluster_atoms,
                format="xyz",
            )

    return quantum_cluster_idx, ecp_region_idx
