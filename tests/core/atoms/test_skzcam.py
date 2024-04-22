from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from quacc.atoms.skzcam import (
    _find_cation_shells,
    _get_anion_coordination,
    _get_atom_distances,
    _get_ecp_region,
    convert_pun_to_atoms,
    create_skzcam_clusters,
    get_cluster_info_from_slab,
    insert_adsorbate_to_embedded_cluster,
)

FILE_DIR = Path(__file__).parent


@pytest.fixture()
def embedded_cluster():
    return convert_pun_to_atoms(
        Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
        {"Mg": 2.0, "O": -2.0},
    )


@pytest.fixture()
def distance_matrix(embedded_cluster):
    return embedded_cluster.get_all_distances()


def test_get_cluster_info_from_slab():
    (
        adsorbate,
        slab,
        slab_first_atom_idx,
        center_position,
        adsorbate_vector_from_slab,
    ) = get_cluster_info_from_slab(
        Path(FILE_DIR, "skzcam_files", "NO_MgO.poscar.gz"),
        adsorbate_idx=[0, 1],
        slab_center_idx=[32, 33],
    )

    # Check adsorbate matches reference
    np.testing.assert_allclose(
        adsorbate.get_positions(),
        np.array(
            [
                [5.39130495, 4.07523845, 15.96981134],
                [5.88635842, 4.84892196, 16.72270959],
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )
    assert np.all(adsorbate.get_atomic_numbers().tolist() == [7, 8])

    # Check slab matches reference
    np.testing.assert_allclose(
        slab.get_positions()[::10],
        np.array(
            [
                [0.0, 6.33073849, 7.5],
                [0.0, 4.22049233, 9.61024616],
                [4.2206809, 6.32743192, 11.73976183],
                [4.2019821, 4.21892378, 13.89202884],
                [0.0, 2.11024616, 9.61024616],
                [4.22049233, 0.0, 7.5],
                [4.22098271, 2.10239745, 13.86181098],
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )
    assert np.all(
        slab.get_atomic_numbers().tolist()
        == [
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            12,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
            8,
        ]
    )

    # Check first atom index of slab
    assert slab_first_atom_idx == 30

    # Check center_position matches reference
    np.testing.assert_allclose(
        center_position,
        np.array([1.06307888, -1.06176564, 2.4591779]),
        rtol=1e-05,
        atol=1e-07,
    )

    # Check vector distance of adsorbate from first center atom (corresponding to first atom index) of slab matches reference
    np.testing.assert_allclose(
        adsorbate_vector_from_slab,
        np.array([1.18932285, -0.14368533, 2.0777825]),
        rtol=1e-05,
        atol=1e-07,
    )


def test_generate_chemshell_cluster():
    from quacc.atoms.skzcam import generate_chemshell_cluster

    # First create the slab
    slab = read(Path(FILE_DIR, "skzcam_files", "NO_MgO.poscar.gz"))
    slab = slab[2:].copy()

    # Run ChemShell
    generate_chemshell_cluster(
        slab,
        30,
        {"Mg": 2.0, "O": -2.0},
        filepath=Path(FILE_DIR, "ChemShell_cluster"),
        chemsh_radius_active=15.0,
        chemsh_radius_cluster=25.0,
        write_xyz_file=True,
    )

    # Read the output .xyz file
    chemshell_embedded_cluster = read(Path(FILE_DIR, "ChemShell_cluster.xyz"))

    # Remove from original folder
    if os.path.isfile(Path(FILE_DIR, "ChemShell_cluster.xyz")):
        os.remove(Path(FILE_DIR, "ChemShell_cluster.xyz"))

    # Check that the positions and atomic numbers match reference
    np.testing.assert_allclose(
        chemshell_embedded_cluster.get_positions()[::100],
        np.array(
            [
                [0.00000000e00, 0.00000000e00, 0.00000000e00],
                [-2.09173593e00, -2.10867761e00, -6.39202884e00],
                [2.12875640e00, -6.32916994e00, -6.39202884e00],
                [-2.09273725e00, 1.05516878e01, -2.16301583e00],
                [2.12875640e00, -1.05496623e01, -6.39202884e00],
                [6.34924872e00, -1.05496623e01, -6.39202884e00],
                [1.05725789e01, -1.05444085e01, -2.15965963e00],
                [1.47875715e01, 6.33408913e00, -2.16464681e00],
                [6.34924872e00, -1.47701546e01, -6.39202884e00],
                [1.69010014e01, 6.33551965e00, -2.15224877e00],
                [1.05697410e01, -1.47701546e01, -6.39202884e00],
                [1.05637735e01, 1.68825241e01, -2.17052139e00],
                [-1.68651820e01, 1.26649992e01, -5.68710477e-02],
                [-1.89763671e01, -1.05478802e01, -2.16464681e00],
                [1.05697410e01, -1.89906469e01, -6.39202884e00],
                [-2.31906127e01, -4.21607826e00, -1.24998430e-02],
                [1.47951600e01, 1.89994594e01, -5.11097275e-02],
                [-2.31941976e01, -6.32916994e00, -6.39202884e00],
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        chemshell_embedded_cluster.get_atomic_numbers()[::20].tolist()
        == [
            12,
            12,
            12,
            8,
            12,
            8,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            8,
            8,
            8,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            12,
            12,
            8,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            8,
            12,
            12,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            8,
            8,
            12,
            8,
            12,
            8,
            12,
            12,
            8,
            12,
            12,
            8,
            8,
            8,
            12,
            12,
            8,
            8,
            8,
            12,
            8,
            12,
            8,
            8,
            12,
            12,
            8,
            12,
            8,
            12,
            12,
            8,
            8,
            12,
            9,
            9,
            9,
        ]
    )


def test_convert_pun_to_atoms():
    embedded_cluster = convert_pun_to_atoms(
        Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
        {"Mg": 2.0, "O": -2.0},
    )

    # Check that number of atoms matches our reference
    assert len(embedded_cluster) == 390

    # Check that last 10 elements of the oxi_state match our reference
    np.testing.assert_allclose(
        embedded_cluster.get_array("oxi_states")[-10:],
        np.array(
            [
                -0.80812511,
                2.14427889,
                -0.96000248,
                2.14427887,
                -0.8081251,
                2.10472993,
                -0.89052904,
                2.10472993,
                -0.8081251,
                2.14427887,
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )

    # Check that first 10 elements of atom_type array match our reference
    assert np.all(
        embedded_cluster.get_array("atom_type")[:10]
        == [
            "cation",
            "anion",
            "anion",
            "anion",
            "anion",
            "anion",
            "cation",
            "cation",
            "cation",
            "cation",
        ]
    )

    # Check that the positions of the atom matches
    np.testing.assert_allclose(
        embedded_cluster[200].position,
        np.array([6.33074029, -2.11024676, -6.37814205]),
        rtol=1e-05,
        atol=1e-07,
    )


def test_insert_adsorbate_to_embedded_cluster(embedded_cluster):
    # Create a CO molecule
    adsorbate = Atoms(
        "CO", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.128]], pbc=[False, False, False]
    )

    # Insert the CO molecule to the embedded cluster
    embedded_cluster, quantum_idx, ecp_idx = insert_adsorbate_to_embedded_cluster(
        embedded_cluster,
        adsorbate,
        [0.0, 0.0, 2.0],
        [[0, 1, 3, 4], [5, 6, 7, 8]],
        [[0, 1, 3, 4], [5, 6, 7, 8]],
    )

    # Check that the positions of the first 10 atoms of the embedded cluster matches the reference positions, oxi_states and atom_type
    np.testing.assert_allclose(
        embedded_cluster.get_positions()[:10],
        np.array(
            [
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 3.128],
                [0.0, 0.0, 0.0],
                [-2.12018426, 0.0, 0.00567209],
                [0.0, 2.12018426, 0.00567209],
                [2.12018426, 0.0, 0.00567209],
                [0.0, -2.12018426, 0.00567209],
                [0.0, 0.0, -2.14129966],
                [-2.11144262, 2.11144262, -0.04367284],
                [2.11144262, 2.11144262, -0.04367284],
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        embedded_cluster.get_chemical_symbols()[:10]
        == ["C", "O", "Mg", "O", "O", "O", "O", "O", "Mg", "Mg"]
    )
    np.testing.assert_allclose(
        embedded_cluster.get_array("oxi_states")[:10],
        np.array([0.0, 0.0, 2.0, -2.0, -2.0, -2.0, -2.0, -2.0, 2.0, 2.0]),
        rtol=1e-05,
        atol=1e-07,
    )
    assert np.all(
        embedded_cluster.get_array("atom_type")[:10]
        == [
            "adsorbate",
            "adsorbate",
            "cation",
            "anion",
            "anion",
            "anion",
            "anion",
            "anion",
            "cation",
            "cation",
        ]
    )

    # Check that the quantum_idx and ecp_idx match the reference
    assert np.all(quantum_idx == [[2, 3, 5, 6], [7, 8, 9, 10]])
    assert np.all(ecp_idx == [[2, 3, 5, 6], [7, 8, 9, 10]])


def test_get_atom_distances():
    # Creating a H2 molecule as an Atoms object
    h2_molecule = Atoms("H2", positions=[(0, 0, 0), (0, 0, 2)])

    # Run _get_atom_distances function to get distance of h2 molecule atoms from a center position
    atom_distances = _get_atom_distances(h2_molecule, [2, 0, 0])

    np.testing.assert_allclose(
        atom_distances, np.array([2.0, 2.82842712]), rtol=1e-05, atol=1e-07
    )


def test_find_cation_shells(embedded_cluster):
    # Get distance of atoms from the center
    distances = _get_atom_distances(embedded_cluster, [0, 0, 2])

    # Find the cation shells from the distances
    cation_shells, cation_shells_idx = _find_cation_shells(
        embedded_cluster, distances, shell_width=0.005
    )

    # As these list of lists do not have the same length, we flatten first 5 lists into a 1D list for comparison
    cation_shells_flatten = [item for row in cation_shells[:5] for item in row]
    cation_shells_idx_flatten = [item for row in cation_shells_idx[:5] for item in row]

    # Check that these lists are correct
    np.testing.assert_allclose(
        cation_shells_flatten,
        np.array(
            [
                2.0,
                3.6184221134101624,
                3.6184221134101655,
                3.6184221134101655,
                3.6184221134101686,
                4.646732760541734,
                4.646732760541734,
                4.646732760541736,
                4.646732760541736,
                4.6888354582307805,
                4.6888354582307805,
                4.6888354582307805,
                4.6888354582307805,
                6.267895285274443,
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        cation_shells_idx_flatten == [0, 9, 8, 6, 7, 11, 12, 10, 13, 19, 21, 18, 20, 22]
    )


def test_get_anion_coordination(embedded_cluster, distance_matrix):
    # Get the anions for the second SKZCAM shell
    anion_shell_idx = _get_anion_coordination(
        embedded_cluster, [9, 8, 6, 7], distance_matrix
    )

    # Check anion indices match with reference
    assert np.all(
        anion_shell_idx == [1, 2, 3, 4, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 29, 30]
    )


def test_get_ecp_region(embedded_cluster, distance_matrix):
    # Find the ECP region for the first cluster
    ecp_region_idx = _get_ecp_region(
        embedded_cluster,
        quantum_cluster_idx=[[0, 1, 2, 3, 4, 5]],
        dist_matrix=distance_matrix,
        ecp_dist=3,
    )

    # Check ECP region indices match with reference
    assert np.all(ecp_region_idx[0] == [6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22])


def test_create_skzcam_clusters(tmpdir):
    # Get quantum cluster and ECP region indices
    embedded_cluster, quantum_cluster_idx, ecp_region_idx = create_skzcam_clusters(
        Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
        [0, 0, 2],
        {"Mg": 2.0, "O": -2.0},
        shell_max=2,
        ecp_dist=3.0,
        shell_width=0.005,
        write_clusters=True,
        write_clusters_path=Path(tmpdir),
    )

    # Check quantum cluster indices match with reference
    assert np.all(
        quantum_cluster_idx[1]
        == [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            14,
            15,
            16,
            17,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
        ]
    )

    # Check ECP region indices match with reference
    assert np.all(
        ecp_region_idx[1]
        == [
            10,
            11,
            12,
            13,
            18,
            19,
            20,
            21,
            22,
            39,
            40,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            54,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            83,
        ]
    )
    # Read the written output and check that it matches with the reference positions and atomic numbers
    skzcam_cluster = read(Path(tmpdir, "SKZCAM_cluster_0.xyz"))

    np.testing.assert_allclose(
        skzcam_cluster.get_positions(),
        np.array(
            [
                [0.0, 0.0, 0.0],
                [-2.12018426, 0.0, 0.00567209],
                [0.0, 2.12018426, 0.00567209],
                [2.12018426, 0.0, 0.00567209],
                [0.0, -2.12018426, 0.00567209],
                [0.0, 0.0, -2.14129966],
            ]
        ),
        rtol=1e-04,
        atol=1e-07,
    )

    assert np.all(skzcam_cluster.get_atomic_numbers().tolist() == [12, 8, 8, 8, 8, 8])
