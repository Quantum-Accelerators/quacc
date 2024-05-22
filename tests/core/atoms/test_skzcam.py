from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io import read
from numpy.testing import assert_allclose

from quacc.atoms.skzcam import (
    _find_cation_shells,
    _get_anion_coordination,
    _get_atom_distances,
    _get_ecp_region,
    convert_pun_to_atoms,
    create_atom_coord_string,
    create_orca_eint_blocks,
    create_orca_point_charge_file,
    create_skzcam_clusters,
    format_ecp_info,
    generate_coords_block,
    generate_orca_input_preamble,
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
def embedded_adsorbed_cluster():
    embedded_cluster, quantum_cluster_indices, ecp_region_indices = (
        create_skzcam_clusters(
            Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
            [0, 0, 2],
            {"Mg": 2.0, "O": -2.0},
            shell_max=2,
            ecp_dist=3,
            write_clusters=False,
        )
    )
    adsorbate = Atoms(
        "CO", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.128]], pbc=[False, False, False]
    )

    embedded_adsorbed_cluster, quantum_cluster_indices, ecp_region_indices = (
        insert_adsorbate_to_embedded_cluster(
            embedded_cluster,
            adsorbate,
            [0.0, 0.0, 2.0],
            quantum_cluster_indices,
            ecp_region_indices,
        )
    )
    return embedded_adsorbed_cluster


@pytest.fixture()
def distance_matrix(embedded_cluster):
    return embedded_cluster.get_all_distances()


def test_create_orca_eint_blocks(embedded_adsorbed_cluster):
    element_info = {
        "C": {
            "basis": "aug-cc-pVDZ",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "aug-cc-pVDZ/C",
        },
        "O": {
            "basis": "aug-cc-pVDZ",
            "core": 2,
            "ri_scf_basis": "def2/JK",
            "ri_cwft_basis": "aug-cc-pVDZ/C",
        },
        "Mg": {
            "basis": "cc-pVDZ",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "cc-pVDZ/C",
        },
    }

    pal_nprocs_block = {"nprocs": 1, "maxcore": 5000}

    method_block = {"Method": "hf", "RI": "on", "RunTyp": "Energy"}

    scf_block = {
        "HFTyp": "rhf",
        "Guess": "MORead",
        "MOInp": '"orca_svp_start.gbw"',
        "SCFMode": "Direct",
        "sthresh": "1e-6",
        "AutoTRAHIter": 60,
        "MaxIter": 1000,
    }

    ecp_info = {
        "Mg": """NewECP
N_core 0
lmax f
s 1
1      1.732000000   14.676000000 2
p 1
1      1.115000000    5.175700000 2
d 1
1      1.203000000   -1.816000000 2
f 1
1      1.000000000    0.000000000 2
end"""
    }

    ad_slab_block, ad_block, slab_block = create_orca_eint_blocks(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info,
        pal_nprocs_block,
        method_block,
        scf_block,
        ecp_info,
        include_cp=True,
        multiplicity={"ad_slab": 1, "slab": 2, "ad": 3},
    )
    # Check that the strings and floats in ad_slab_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in ad_slab_block.split()[::10]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [
            3.128,
            0.0,
            0.00567209089,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
        ],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in ad_slab_block.split()[::5]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%pal",
            "Method",
            "Energy",
            "NewNCore",
            "O",
            "NewGTO",
            "Mg",
            '"aug-cc-pVDZ"',
            "end",
            "NewAuxJGTO",
            "C",
            '"cc-pVDZ/C"',
            "end",
            "Guess",
            "Direct",
            "MaxIter",
            "xyz",
            "Charge",
            "O",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "end",
        ]
    )

    # Check that the strings and floats in ad_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in ad_block.split()[::2]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [
            1.0,
            2.0,
            2.0,
            2.0,
            0.0,
            2.0,
            0.0,
            3.128,
            0.0,
            0.0,
            -2.12018425659,
            0.00567209089,
            0.0,
            0.00567209089,
            2.12018425659,
            0.00567209089,
            0.0,
            0.00567209089,
            0.0,
            -2.14129966123,
        ],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in ad_block.split()[::2]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%pal",
            "%maxcore",
            "end",
            '"orca.pc"',
            "Method",
            "RI",
            "RunTyp",
            "NewNCore",
            "NewNCore",
            "NewNCore",
            "end",
            "NewGTO",
            '"aug-cc-pVDZ"',
            "NewGTO",
            '"cc-pVDZ"',
            "NewGTO",
            '"aug-cc-pVDZ"',
            "NewAuxJGTO",
            '"def2/J"',
            "NewAuxJGTO",
            '"def2/J"',
            "NewAuxJGTO",
            '"def2/JK"',
            "NewAuxCGTO",
            '"aug-cc-pVDZ/C"',
            "NewAuxCGTO",
            '"cc-pVDZ/C"',
            "NewAuxCGTO",
            '"aug-cc-pVDZ/C"',
            "end",
            "HFTyp",
            "Guess",
            "MOInp",
            "SCFMode",
            "sthresh",
            "AutoTRAHIter",
            "MaxIter",
            "end",
            "CTyp",
            "Mult",
            "Units",
            "Charge",
            "coords",
            "end",
        ]
    )

    # Check that the strings and floats in slab_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in slab_block.split()[::10]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [
            3.128,
            0.0,
            0.00567209089,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
            5.1757,
            1.0,
            2.0,
            1.203,
        ],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in slab_block.split()[::5]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%pal",
            "Method",
            "Energy",
            "NewNCore",
            "O",
            "NewGTO",
            "Mg",
            '"aug-cc-pVDZ"',
            "end",
            "NewAuxJGTO",
            "C",
            '"cc-pVDZ/C"',
            "end",
            "Guess",
            "Direct",
            "MaxIter",
            "xyz",
            "Charge",
            "O",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "Mg>",
            "NewECP",
            "s",
            "end",
        ]
    )


def test_create_atom_coord_string(embedded_adsorbed_cluster):
    atom = embedded_adsorbed_cluster[0]

    # First let's try the case where it's a normal atom.
    atom_coord_string = create_atom_coord_string(atom)

    with pytest.raises(ValueError):
        create_atom_coord_string(
            atom, atom_ecp_info="NewECP\nECP_info1\nECP_info2\n", ghost_atom=True
        )

    with pytest.raises(ValueError):
        create_atom_coord_string(atom, atom_ecp_info="NewECP\nECP_info1\nECP_info2\n")

    assert (
        atom_coord_string
        == "C                       0.00000000000    0.00000000000    2.00000000000\n"
    )

    # Let's now try the case where it is a ghost atom.
    atom_coord_string = create_atom_coord_string(atom, ghost_atom=True)
    assert (
        atom_coord_string
        == "C:                      0.00000000000    0.00000000000    2.00000000000\n"
    )

    # Let's now try the case where it is an atom in the ECP region.
    atom_coord_string = create_atom_coord_string(
        atom, atom_ecp_info="NewECP\nECP_info1\nECP_info2\n", pc_charge=2.0
    )
    assert (
        atom_coord_string
        == "C>     2.00000000000    0.00000000000    0.00000000000    2.00000000000\nNewECP\nECP_info1\nECP_info2\n"
    )


def test_get_coords_block(embedded_adsorbed_cluster):
    ecp_info = {
        "Mg": """NewECP
N_core 0
lmax f
s 1
1      1.732000000   14.676000000 2
p 1
1      1.115000000    5.175700000 2
d 1
1      1.203000000   -1.816000000 2
f 1
1      1.000000000    0.000000000 2
end"""
    }

    # Confirm that multiplicity given as 1 if not provided
    ad_slab_coords, ad_coords, slab_coords = generate_coords_block(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        ecp_info,
    )

    assert ad_slab_coords.split()[4] == "1"
    assert slab_coords.split()[4] == "1"
    assert ad_coords.split()[4] == "1"

    ad_slab_coords, ad_coords, slab_coords = generate_coords_block(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        ecp_info,
        include_cp=True,
        multiplicity={"ad_slab": 1, "slab": 2, "ad": 3},
    )

    # Check that the strings and floats in ad_slab_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in ad_slab_coords.split()[::10]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [
            0.0,
            0.0,
            1.0,
            2.11144262254,
            1.732,
            1.0,
            2.0,
            1.0,
            -2.11144262254,
            1.732,
            1.0,
            2.0,
            1.0,
            2.10705287155,
            1.732,
            1.0,
            2.0,
            1.0,
            -2.10705287155,
            1.732,
            1.0,
            2.0,
            1.0,
            4.22049352791,
            1.732,
            1.0,
            2.0,
            1.0,
            -4.22049352791,
            1.732,
            1.0,
            2.0,
            1.0,
        ],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in ad_slab_coords.split()[::5]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%coords",
            "Units",
            "C",
            "O",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
        ]
    )

    # Check that the strings and floats in ad_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in ad_coords.split()[::2]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.12018425659, 0.0, -2.12018425659, 0.0],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in ad_coords.split()[::2]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%coords",
            "xyz",
            "angs",
            "C",
            "O",
            "Mg:",
            "O:",
            "O:",
            "O:",
            "O:",
            "O:",
            "end",
        ]
    )

    # Check that the strings and floats in slab_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in slab_coords.split()[::10]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [
            0.0,
            0.0,
            1.0,
            2.11144262254,
            1.732,
            1.0,
            2.0,
            1.0,
            -2.11144262254,
            1.732,
            1.0,
            2.0,
            1.0,
            2.10705287155,
            1.732,
            1.0,
            2.0,
            1.0,
            -2.10705287155,
            1.732,
            1.0,
            2.0,
            1.0,
            4.22049352791,
            1.732,
            1.0,
            2.0,
            1.0,
            -4.22049352791,
            1.732,
            1.0,
            2.0,
            1.0,
        ],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in slab_coords.split()[::5]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%coords",
            "Units",
            "C:",
            "O",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
            "lmax",
            "f",
        ]
    )

    # Also check the case where include_cp is False
    ad_slab_coords, ad_coords, slab_coords = generate_coords_block(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        ecp_info,
        include_cp=False,
        multiplicity={"ad_slab": 1, "slab": 2, "ad": 3},
    )

    # Check that the strings and floats in ad_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in ad_coords.split()[::2]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [3.0, 0.0, 0.0, 0.0],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in ad_coords.split()[::2]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == ["%coords", "xyz", "angs", "C", "O", "end"]
    )

    # Check that the strings and float in slab_coords matches reference
    assert_allclose(
        [
            float(x)
            for x in slab_coords.split()[::10]
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ],
        [
            2.12018425659,
            -1.816,
            2.0,
            1.0,
            2.0,
            1.0,
            -1.816,
            2.0,
            1.0,
            2.0,
            1.0,
            -1.816,
            2.0,
            1.0,
            2.0,
            1.0,
            -1.816,
            2.0,
            1.0,
            2.0,
            1.0,
            -1.816,
            2.0,
            1.0,
            2.0,
            1.0,
            -1.816,
            2.0,
            1.0,
            2.0,
            1.0,
            -1.816,
        ],
        rtol=1e-05,
        atol=1e-07,
    )

    assert np.all(
        [
            x
            for x in slab_coords.split()[::5]
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ]
        == [
            "%coords",
            "Units",
            "Mg",
            "O",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "N_core",
            "p",
            "end",
        ]
    )


def test_format_ecp_info():
    atom_ecp_info = """dummy_info

NewECP


N_core 0
lmax s
s 1
1      1.732000000   14.676000000 2
eNd

"""
    formatted_atom_ecp_info = format_ecp_info(atom_ecp_info)

    with pytest.raises(ValueError):
        format_ecp_info("dummy_info\nN_core0\nend")

    assert (
        formatted_atom_ecp_info
        == "NewECP\nN_core 0\nlmax s\ns 1\n1      1.732000000   14.676000000 2\nend\n"
    )

    atom_ecp_info = """
N_core 0
lmax s
s 1
1      1.732000000   14.676000000 2

"""
    assert (
        formatted_atom_ecp_info
        == "NewECP\nN_core 0\nlmax s\ns 1\n1      1.732000000   14.676000000 2\nend\n"
    )


def test_generate_orca_input_preamble(embedded_adsorbed_cluster):
    # Set-up some information needed for generating orca input
    element_info = {
        "C": {
            "basis": "aug-cc-pVDZ",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "aug-cc-pVDZ/C",
        },
        "O": {
            "basis": "aug-cc-pVDZ",
            "core": 2,
            "ri_scf_basis": "def2/JK",
            "ri_cwft_basis": "aug-cc-pVDZ/C",
        },
        "Mg": {
            "basis": "cc-pVDZ",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "cc-pVDZ/C",
        },
    }

    pal_nprocs_block = {"nprocs": 1, "maxcore": 5000}

    method_block = {"Method": "hf", "RI": "on", "RunTyp": "Energy"}

    scf_block = {
        "HFTyp": "rhf",
        "Guess": "MORead",
        "MOInp": '"orca_svp_start.gbw"',
        "SCFMode": "Direct",
        "sthresh": "1e-6",
        "AutoTRAHIter": 60,
        "MaxIter": 1000,
    }

    # Check whether error raised if not all element_info is not provided
    with pytest.raises(ValueError):
        element_info_error = {"C": element_info["C"]}
        generate_orca_input_preamble(
            embedded_adsorbed_cluster,
            [0, 1, 2, 3, 4, 5, 6, 7],
            element_info=element_info_error,
            pal_nprocs_block=pal_nprocs_block,
            method_block=method_block,
            scf_block=scf_block,
        )

    # Generate the orca input preamble
    preamble_input = generate_orca_input_preamble(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        element_info=element_info,
        pal_nprocs_block=pal_nprocs_block,
        method_block=method_block,
        scf_block=scf_block,
    )

    assert (
        preamble_input
        == '%pal nprocs 1 end\n%maxcore 5000 end\n%pointcharges "orca.pc"\n%method\nMethod hf\nRI on\nRunTyp Energy\nNewNCore C 2 end\nNewNCore Mg 2 end\nNewNCore O 2 end\nend\n%basis\nNewGTO C "aug-cc-pVDZ" end\nNewGTO Mg "cc-pVDZ" end\nNewGTO O "aug-cc-pVDZ" end\nNewAuxJGTO C "def2/J" end\nNewAuxJGTO Mg "def2/J" end\nNewAuxJGTO O "def2/JK" end\nNewAuxCGTO C "aug-cc-pVDZ/C" end\nNewAuxCGTO Mg "cc-pVDZ/C" end\nNewAuxCGTO O "aug-cc-pVDZ/C" end\nend\n%scf\nHFTyp rhf\nGuess MORead\nMOInp "orca_svp_start.gbw"\nSCFMode Direct\nsthresh 1e-6\nAutoTRAHIter 60\nMaxIter 1000\nend\n'
    )

    # Check the case if the element_info has all of the same values
    element_info = {
        "C": {
            "basis": "def2-SVP",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "def2-SVP/C",
        },
        "O": {
            "basis": "def2-SVP",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "def2-SVP/C",
        },
        "Mg": {
            "basis": "def2-SVP",
            "core": 2,
            "ri_scf_basis": "def2/J",
            "ri_cwft_basis": "def2-SVP/C",
        },
    }

    preamble_input = generate_orca_input_preamble(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        element_info=element_info,
        pal_nprocs_block=pal_nprocs_block,
        method_block=method_block,
        scf_block=scf_block,
    )

    assert (
        preamble_input
        == '%pal nprocs 1 end\n%maxcore 5000 end\n%pointcharges "orca.pc"\n%method\nMethod hf\nRI on\nRunTyp Energy\nNewNCore C 2 end\nNewNCore Mg 2 end\nNewNCore O 2 end\nend\n%basis\nBasis def2-SVP\nAux def2/J\nAuxC def2-SVP/C\nend\n%scf\nHFTyp rhf\nGuess MORead\nMOInp "orca_svp_start.gbw"\nSCFMode Direct\nsthresh 1e-6\nAutoTRAHIter 60\nMaxIter 1000\nend\n'
    )

    # Testing the case if we provide no blocks
    preamble_input = generate_orca_input_preamble(
        embedded_adsorbed_cluster, [0, 1, 2, 3, 4, 5, 6, 7]
    )

    assert preamble_input == '%pointcharges "orca.pc"\n'


def test_create_orca_point_charge_file(embedded_adsorbed_cluster, tmpdir):
    # Create the point charge file
    create_orca_point_charge_file(
        embedded_adsorbed_cluster,
        [0, 1, 2, 3, 4, 5, 6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        Path(tmpdir, "orca.pc"),
    )

    # Test whether exception is raised if indices shared between quantum region and ecp region
    with pytest.raises(ValueError):
        create_orca_point_charge_file(
            embedded_adsorbed_cluster,
            [0, 1, 2, 3, 4, 5, 6, 7],
            [7, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
            Path(tmpdir, "orca.pc"),
        )

    # Read the written file
    orca_pc_file = np.loadtxt(Path(tmpdir, "orca.pc"), skiprows=1)

    # Check that the contents of the file match the reference
    assert len(orca_pc_file) == 371

    assert_allclose(
        orca_pc_file[::30],
        np.array(
            [
                [-2.00000000e00, -2.11070451e00, 2.11070451e00, -2.14923990e00],
                [2.00000000e00, 2.11024676e00, -2.11024676e00, -4.26789529e00],
                [2.00000000e00, 6.32954443e00, 2.11144262e00, -4.36728442e-02],
                [-2.00000000e00, -4.22049353e00, 6.32889566e00, 7.72802266e-03],
                [2.00000000e00, -6.33074029e00, -2.11024676e00, -4.26789529e00],
                [-2.00000000e00, 4.22049353e00, -6.33074029e00, -4.26789529e00],
                [-2.00000000e00, 6.33074029e00, 2.11024676e00, -6.37814205e00],
                [-2.00000000e00, 2.11024676e00, -8.44098706e00, -4.26789529e00],
                [-2.00000000e00, -8.44098706e00, -6.32080280e00, 5.67209089e-03],
                [2.00000000e00, -2.11024676e00, 8.44098706e00, -6.37814205e00],
                [8.00000000e-01, -4.64254288e01, 3.79844418e01, -3.99237095e-02],
                [3.12302613e00, -0.00000000e00, -5.71441194e01, -2.36698692e01],
                [2.10472999e00, -2.36698692e01, 5.71441194e01, 2.59086514e01],
            ]
        ),
        rtol=1e-05,
        atol=1e-07,
    )


def test_get_cluster_info_from_slab():
    (
        adsorbate,
        slab,
        slab_first_atom_idx,
        center_position,
        adsorbate_vector_from_slab,
    ) = get_cluster_info_from_slab(
        Path(FILE_DIR, "skzcam_files", "NO_MgO.poscar.gz"), [32, 33], [0, 1]
    )

    # Check adsorbate matches reference
    assert_allclose(
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
    assert_allclose(
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
    assert_allclose(
        center_position,
        np.array([1.06307888, -1.06176564, 2.4591779]),
        rtol=1e-05,
        atol=1e-07,
    )

    # Check vector distance of adsorbate from first center atom (corresponding to first atom index) of slab matches reference
    assert_allclose(
        adsorbate_vector_from_slab,
        np.array([1.18932285, -0.14368533, 2.0777825]),
        rtol=1e-05,
        atol=1e-07,
    )


def test_generate_chemshell_cluster(tmp_path):
    from quacc.atoms.skzcam import generate_chemshell_cluster

    # First create the slab
    slab = read(Path(FILE_DIR, "skzcam_files", "NO_MgO.poscar.gz"))[2:]

    # Run ChemShell
    generate_chemshell_cluster(
        slab,
        30,
        {"Mg": 2.0, "O": -2.0},
        tmp_path,
        chemsh_radius_active=15.0,
        chemsh_radius_cluster=25.0,
        write_xyz_file=True,
    )

    # Read the output .xyz file
    chemshell_embedded_cluster = read(tmp_path / "ChemShell_cluster.xyz")

    # Check that the positions and atomic numbers match reference
    assert_allclose(
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
    assert_allclose(
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
    assert_allclose(
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
    assert_allclose(
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
    assert_allclose(
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
    assert np.all(quantum_idx == [[0, 1, 2, 3, 5, 6], [0, 1, 7, 8, 9, 10]])
    assert np.all(ecp_idx == [[2, 3, 5, 6], [7, 8, 9, 10]])


def test_get_atom_distances():
    # Creating a H2 molecule as an Atoms object
    h2_molecule = Atoms("H2", positions=[(0, 0, 0), (0, 0, 2)])

    # Run _get_atom_distances function to get distance of h2 molecule atoms from a center position
    atom_distances = _get_atom_distances(h2_molecule, [2, 0, 0])

    assert_allclose(atom_distances, np.array([2.0, 2.82842712]), rtol=1e-05, atol=1e-07)


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
    assert_allclose(
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
        embedded_cluster, [[0, 1, 2, 3, 4, 5]], distance_matrix, ecp_dist=3
    )

    # Check ECP region indices match with reference
    assert np.all(ecp_region_idx[0] == [6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22])


def test_create_skzcam_clusters(tmp_path):
    # Get quantum cluster and ECP region indices
    _, quantum_cluster_idx, ecp_region_idx = create_skzcam_clusters(
        Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
        [0, 0, 2],
        {"Mg": 2.0, "O": -2.0},
        shell_max=2,
        ecp_dist=3.0,
        shell_width=0.005,
        write_clusters=True,
        write_clusters_path=tmp_path,
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
    skzcam_cluster = read(tmp_path / "SKZCAM_cluster_0.xyz")

    assert_allclose(
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
