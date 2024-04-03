from __future__ import annotations

from pathlib import Path
from typing import Dict, Union

import numpy as np
from ase import Atoms, io
from ase.units import Bohr


def convert_pun_to_Atoms(
    pun_filename: Union[str, Path], atom_oxi_states: Dict
) -> Atoms:
    """
    Reads a .pun file and returns an ASE Atoms object containing the atomic coordinates, point charges/oxidation states and atom types.

    Parameters
    ----------
    pun_filename : str
        The path to the .pun file created by ChemShell to be read.
    atom_oxi_states : dict
        A dictionary containing the atomic symbols as keys and the oxidation states as values.

    Returns
    -------
    Atoms
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file.
    """

    # Create an emtpy emb_clus Atoms object
    emb_clus = Atoms()

    # Create a dictionary containing the atom types and whether they are cations or anions
    atom_type_dict = {}
    for atom, oxi_state in atom_oxi_states.items():
        if oxi_state > 0:
            atom_type_dict[atom] = "cation"
        elif oxi_state < 0:
            atom_type_dict[atom] = "anion"
        else:
            atom_type_dict[atom] = "neutral"

    # Load the pun file as a list of strings
    raw_pun_file = [line.rstrip() for line in open(pun_filename)]

    # Get the number of atoms and number of atomic charges in the .pun file
    num_atoms = int(raw_pun_file[3].split()[-1])
    num_charges = int(raw_pun_file[4 + num_atoms - 1 + 3].split()[-1])

    # Check if number of atom charges same as number of atom positions
    if num_atoms != num_charges:
        raise ValueError(
            "Number of atomic positions and atomic charges in the .pun file are not the same."
        )

    raw_atom_positions = raw_pun_file[4 : 4 + num_atoms]
    raw_charge_list = raw_pun_file[7 + num_atoms : 7 + 2 * num_atoms]
    charge_list = [float(charge) for charge in raw_charge_list]
    atom_type_list = []
    # Add the atomic positions the emb_clus Atoms object (converting from Bohr to Angstrom)
    for i, line in enumerate(raw_atom_positions):
        line_info = line.split()
        # Add the atom type to the atom_type_list
        if line_info[0] not in atom_type_dict:
            if line_info[0] == "F":
                atom_type_list.append("fitting charge")
            else:
                atom_type_list.append("unknown")
        else:
            atom_type_list.append(atom_type_dict[line_info[0]])

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
        emb_clus = emb_clus.copy() + atom.copy()

    # Centre the embedded cluster so that atom 1 is at the [0, 0, 0] position
    emb_clus.translate(-emb_clus[0].position)

    # Add the charge_list to the Atoms object
    emb_clus.set_array("oxi_states", np.array(charge_list))
    emb_clus.set_array("atom_type", np.array(atom_type_list))

    return emb_clus


def get_atom_distances(emb_clus: Atoms, centre_position: np.ndarray) -> np.ndarray:
    """
    Returns the distance of all atoms from the centre position of the embedded cluster

    Parameters
    ----------
    emb_clus : Atoms
        The ASE Atoms object containing the atomic coordinates.

    Returns
    -------
    np.ndarray
        An array containing the distances of each atom in the Atoms object fro the cluster centre.
    """

    # Get the number of atoms in the Atoms object
    num_atoms = len(emb_clus)

    # Create an empty matrix to store the distances between each pair of atoms
    distances = np.zeros(num_atoms)

    # Calculate the distances between each pair of atoms
    for i in range(num_atoms):
        distances[i] = np.linalg.norm(emb_clus[i].position - centre_position)

    return distances


def find_cation_shells(
    emb_clus: Atoms,
    distances: np.ndarray,
    shell_width: float = 0.005,
) -> list:
    """
    Returns a list of lists containing the indices of the cations in each shell, based on distance from the embedded cluster centre. This is achieved by clustering the data based on the DBSCAN clustering algorithm.

    Parameters
    ----------
    emb_clus : Atoms
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

    distances_sorted = []
    distances_sorted_idx = []
    # Sort the points by distance from the cluster centre for the cations only
    for i in np.argsort(distances):
        if emb_clus.get_array("atom_type")[i] == "cation":
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
    emb_clus: Atoms,
    cation_shell_idx: list,
    dist_matrix: np.ndarray,
    bond_dist: float = 2.5,
) -> list:
    """
    Returns a list of lists containing the indices of the anions coordinating the cation indices provided.

    Parameters
    ----------
    emb_clus : Atoms
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    cation_shell_idx : list
        A list of the indices of the cations in the cluster.
    dist_matrix : np.ndarray
        A matrix containing the distances between each pair of atoms in the embedded cluster.

    Returns
    -------
    anion_coord_idx: list
        A list containing the indices of the anions coordinating the cation indices.
    """

    # Define the empty list to store the anion coordination
    anion_coord_idx = []

    # Iterate over the cation shell indices and find the atoms within the bond distance of each cation
    for atom_idx in cation_shell_idx:
        anion_coord_idx += [
            idx
            for idx, dist in enumerate(dist_matrix[atom_idx])
            if (dist < bond_dist and emb_clus.get_array("atom_type")[idx] == "anion")
        ]

    return list(set(anion_coord_idx))


def get_ecp_region(
    emb_clus: Atoms,
    quantum_cluster_idx: list,
    dist_matrix: np.ndarray,
    ecp_dist: float = 6.0,
) -> list:
    """
    Returns a list of lists containing the indices of the atoms in the ECP region of the embedded cluster for each quantum cluster

    Parameters
    ----------
    emb_clus : Atoms
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    quantum_cluster_idx : list
        A list of lists containing the indices of the atoms in each quantum cluster.
    dist_matrix : np.ndarray
        A matrix containing the distances between each pair of atoms in the embedded cluster.
    ecp_dist : float
        The distance from edges of the quantum cluster to define the ECP region.

    Returns
    -------
    ecp_region_idx: list
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
                    and emb_clus.get_array("atom_type")[idx] == "cation"
                ):
                    cluster_ecp_region_idx += [idx]

        ecp_region_idx += [list(set(cluster_ecp_region_idx))]

    return ecp_region_idx


def create_skzcam_clusters(
    pun_filename: Union[str, Path],
    centre_position: np.ndarray,
    atom_oxi_states: Dict,
    shell_max: int = 10,
    shell_width: float = 0.005,
    bond_dist: float = 2.5,
    write_clusters: bool = False,
    write_clusters_path: Union[str, Path] = None,
    write_cluster_name: str = "SKZCAM_cluster",
    write_include_ecp: bool = False,
) -> list:
    """
    Returns a list of list containing the indices of the atoms (in emb_clus) which form the quantum clusters in the SKZCAM protocol. The number of clusters created is controlled by the rdf_max parameter.

    Parameters
    ----------
    emb_clus : Atoms
        The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
    shell_max : int
        The maximum number of quantum clusters to be created.
    write_clusters : bool
        If True, the quantum clusters will be written to a file.
    write_clusters_path : str
        The path to the file where the quantum clusters will be written.

    Returns
    -------
    clusters: list
        A list of lists containing the indices of the atoms in each quantum cluster.
    """

    # Read the .pun file and create the emb_clus Atoms object
    emb_clus = convert_pun_to_Atoms(pun_filename, atom_oxi_states)

    # Get distances of all atoms from the cluster centre
    atom_centre_distances = get_atom_distances(emb_clus, centre_position)

    # Determine the cation shells from the centre of the embedded cluster
    cation_shells, cation_shells_idx = find_cation_shells(
        emb_clus, atom_centre_distances, shell_width
    )

    # Find the maximum idx of cation shells up to shell_max
    # if shell_max + 7 > len(cation_shells_idx):
    #     shell_max1 = len(cation_shells_idx)
    # else:
    #     shell_max1 = shell_max + 7

    # flatten_idx = [x for xs in cation_shells_idx[:shell_max1] for x in xs]
    # max_idx = max(flatten_idx)

    # Make a smaller embedded cluster with atoms up to max_idx closer to the centre
    # emb_clus_small = emb_clus.copy()[:max_idx]

    # Create the distance matrix for the embedded cluster
    emb_clus_all_dist = emb_clus.get_all_distances()

    anion_coord_idx = []
    # Create the anion coordination list for each cation shell
    for shell_idx in range(shell_max):
        cation_shell = cation_shells_idx[shell_idx]
        anion_coord_idx += [
            get_anion_coordination(emb_clus, cation_shell, emb_clus_all_dist, bond_dist)
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
    ecp_region_idx = get_ecp_region(emb_clus, quantum_cluster_idx, emb_clus_all_dist)

    if write_clusters:
        if write_clusters_path is None:
            print(
                "No path provided to write the quantum clusters. Writing to the current working directory."
            )
            write_clusters_path = Path.cwd()
        # Write the quantum clusters to a file
        for idx, cluster in enumerate(quantum_cluster_idx):
            cluster_atoms = emb_clus[cluster]
            if write_include_ecp:
                ecp_atoms = emb_clus[ecp_region_idx[idx]]
                ecp_atoms.set_chemical_symbols(np.array(["U"] * len(ecp_atoms)))
                cluster_atoms = cluster_atoms.copy() + ecp_atoms.copy()

            io.write(
                Path(f"{write_clusters_path}/{write_cluster_name}_{idx}.xyz"),
                cluster_atoms,
                format="xyz",
            )
            # cluster_atoms.write(write_clusters_path + f"_{write_cluster_name}_{idx}.xyz")

    return quantum_cluster_idx, ecp_region_idx
