from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import compare_atoms
from ase.io import read
from numpy.testing import assert_allclose, assert_equal

from quacc.atoms.skzcam import (
    MRCCInputGenerator,
    ORCAInputGenerator,
    _find_cation_shells,
    _get_anion_coordination,
    _get_atom_distances,
    _get_ecp_region,
    convert_pun_to_atoms,
    create_atom_coord_string,
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
def mrcc_input_generator(embedded_adsorbed_cluster, element_info):
    return MRCCInputGenerator(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        ecp_region_indices=[8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info=element_info,
        include_cp=True,
        multiplicities={"adsorbate_slab": 3, "adsorbate": 1, "slab": 2},
    )


@pytest.fixture()
def orca_input_generator(embedded_adsorbed_cluster, element_info):
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
    return ORCAInputGenerator(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        ecp_region_indices=[8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info=element_info,
        include_cp=True,
        multiplicities={"adsorbate_slab": 3, "adsorbate": 1, "slab": 2},
        pal_nprocs_block=pal_nprocs_block,
        method_block=method_block,
        scf_block=scf_block,
        ecp_info=ecp_info,
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
            embedded_cluster=embedded_cluster,
            adsorbate=adsorbate,
            adsorbate_vector_from_slab=[0.0, 0.0, 2.0],
            quantum_cluster_indices=quantum_cluster_indices,
            ecp_region_indices=ecp_region_indices,
        )
    )
    return embedded_adsorbed_cluster


@pytest.fixture()
def element_info():
    return {
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


@pytest.fixture()
def distance_matrix(embedded_cluster):
    return embedded_cluster.get_all_distances()


def test_MRCCInputGenerator_init(embedded_adsorbed_cluster, element_info):
    # Check what happens if multiplicities is not provided
    mrcc_input_generator = MRCCInputGenerator(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        ecp_region_indices=[8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info=element_info,
        include_cp=True,
    )

    assert mrcc_input_generator.multiplicities == {
        "adsorbate_slab": 1,
        "adsorbate": 1,
        "slab": 1,
    }

    mrcc_input_generator = MRCCInputGenerator(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        ecp_region_indices=[8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info=element_info,
        include_cp=True,
        multiplicities={"adsorbate_slab": 3, "adsorbate": 1, "slab": 2},
    )

    assert not compare_atoms(
        mrcc_input_generator.embedded_adsorbed_cluster, embedded_adsorbed_cluster
    )
    assert_equal(mrcc_input_generator.quantum_cluster_indices, [0, 1, 2, 3, 4, 5, 6, 7])
    assert_equal(mrcc_input_generator.adsorbate_indices, [0, 1])
    assert_equal(mrcc_input_generator.slab_indices, [2, 3, 4, 5, 6, 7])
    assert_equal(
        mrcc_input_generator.ecp_region_indices,
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
    )
    assert mrcc_input_generator.element_info == element_info
    assert mrcc_input_generator.include_cp is True
    assert mrcc_input_generator.multiplicities == {
        "adsorbate_slab": 3,
        "adsorbate": 1,
        "slab": 2,
    }

    # Check if error raise if quantum_cluster_indices and ecp_region_indices overlap

    with pytest.raises(
        ValueError, match="An atom in the quantum cluster is also in the ECP region."
    ):
        mrcc_input_generator = MRCCInputGenerator(
            embedded_adsorbed_cluster=embedded_adsorbed_cluster,
            quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
            ecp_region_indices=[7, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
            element_info=element_info,
            include_cp=True,
        )


def test_MRCCInputGenerator_create_eint_blocks(mrcc_input_generator):
    mrcc_input_generator_nocp = deepcopy(mrcc_input_generator)

    mrcc_input_generator_nocp.include_cp = False
    mrcc_input_generator_nocp.generate_input()

    mrcc_input_generator.generate_input()

    reference_block_collated = {
        "adsorbate_slab": {
            "float": [21.0, -2.0, 2.0, 2.0, 2.0, 0.1474277671],
            "string": ["basis_sm=atomtype", "def2/JK", "O"],
        },
        "adsorbate": {"float": [8.0], "string": ["basis_sm=atomtype", "3,4,5,6,7,8"]},
        "slab": {
            "float": [21.0, -2.0, 2.0, 2.0, 2.0, 0.1474277671],
            "string": ["basis_sm=atomtype", "def2/JK", "O"],
        },
    }

    reference_block_nocp_collated = {
        "adsorbate_slab": {
            "float": [21.0, -2.0, 2.0, 2.0, 2.0, 0.1474277671],
            "string": ["basis_sm=atomtype", "def2/JK", "O"],
        },
        "adsorbate": {"float": [2.0], "string": ["basis_sm=atomtype"]},
        "slab": {
            "float": [
                19.0,
                -4.22049352791,
                4.22049352791,
                4.22049352791,
                2.11024676395,
                -0.0,
            ],
            "string": ["basis_sm=atomtype", "no-basis-set", "Mg"],
        },
    }

    generated_block_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }
    generated_block_nocp_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }

    for system in ["adsorbate_slab", "adsorbate", "slab"]:
        generated_block_collated[system]["float"] = [
            float(x)
            for x in mrcc_input_generator.mrccblocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::300]
        generated_block_collated[system]["string"] = [
            x
            for x in mrcc_input_generator.mrccblocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::50]

        generated_block_nocp_collated[system]["float"] = [
            float(x)
            for x in mrcc_input_generator_nocp.mrccblocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::300]
        generated_block_nocp_collated[system]["string"] = [
            x
            for x in mrcc_input_generator_nocp.mrccblocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::50]

        assert_equal(
            generated_block_collated[system]["string"],
            reference_block_collated[system]["string"],
        )
        assert_allclose(
            generated_block_collated[system]["float"],
            reference_block_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )

        assert_equal(
            generated_block_nocp_collated[system]["string"],
            reference_block_nocp_collated[system]["string"],
        )
        assert_allclose(
            generated_block_nocp_collated[system]["float"],
            reference_block_nocp_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )


def test_MRCCInputGenerator_generate_basis_ecp_block(mrcc_input_generator):
    mrcc_input_generator_nocp = deepcopy(mrcc_input_generator)

    mrcc_input_generator_nocp.include_cp = False
    mrcc_input_generator_nocp.generate_basis_ecp_block()

    mrcc_input_generator.generate_basis_ecp_block()

    reference_mrcc_blocks_collated = {
        "adsorbate_slab": [
            "basis_sm=atomtype",
            "no-basis-set",
            "no-basis-set",
            "aug-cc-pVDZ",
            "no-basis-set",
            "def2/JK",
            "no-basis-set",
            "aug-cc-pVDZ/C",
            "no-basis-set",
        ],
        "slab": [
            "basis_sm=atomtype",
            "no-basis-set",
            "no-basis-set",
            "aug-cc-pVDZ",
            "no-basis-set",
            "def2/JK",
            "no-basis-set",
            "aug-cc-pVDZ/C",
            "no-basis-set",
        ],
        "adsorbate": ["basis_sm=atomtype", "aug-cc-pVDZ", "def2/JK", "cc-pVDZ/C"],
    }

    reference_mrcc_blocks_nocp_collated = {
        "adsorbate_slab": [
            "basis_sm=atomtype",
            "no-basis-set",
            "no-basis-set",
            "aug-cc-pVDZ",
            "no-basis-set",
            "def2/JK",
            "no-basis-set",
            "aug-cc-pVDZ/C",
            "no-basis-set",
        ],
        "slab": [
            "basis_sm=atomtype",
            "no-basis-set",
            "basis=atomtype",
            "no-basis-set",
            "dfbasis_scf=atomtype",
            "no-basis-set",
            "dfbasis_cor=atomtype",
            "no-basis-set",
        ],
        "adsorbate": ["basis_sm=atomtype", "aug-cc-pVDZ/C"],
    }

    generated_mrcc_blocks_nocp_collated = {
        system: [] for system in ["adsorbate_slab", "slab", "adsorbate"]
    }
    generated_mrcc_blocks_collated = {
        system: [] for system in ["adsorbate_slab", "slab", "adsorbate"]
    }
    for system in ["adsorbate_slab", "adsorbate", "slab"]:
        generated_mrcc_blocks_collated[system] = mrcc_input_generator.mrccblocks[
            system
        ].split()[::10]
        generated_mrcc_blocks_nocp_collated[system] = (
            mrcc_input_generator_nocp.mrccblocks[system].split()[::10]
        )

        assert_equal(
            generated_mrcc_blocks_collated[system],
            reference_mrcc_blocks_collated[system],
        )
        assert_equal(
            generated_mrcc_blocks_nocp_collated[system],
            reference_mrcc_blocks_nocp_collated[system],
        )


def test_MRCCInputGenerator_create_atomtype_basis(mrcc_input_generator):
    generated_basis_block_without_ecp = mrcc_input_generator.create_atomtype_basis(
        quantum_region=mrcc_input_generator.adsorbate_slab_cluster,
        element_basis_info={
            element: mrcc_input_generator.element_info[element]["ri_cwft_basis"]
            for element in mrcc_input_generator.element_info
        },
    )
    generated_basis_block_with_ecp = mrcc_input_generator.create_atomtype_basis(
        quantum_region=mrcc_input_generator.adsorbate_slab_cluster,
        element_basis_info={
            element: mrcc_input_generator.element_info[element]["ri_cwft_basis"]
            for element in mrcc_input_generator.element_info
        },
        ecp_region=mrcc_input_generator.ecp_region,
    )

    reference_basis_block_without_ecp = "aug-cc-pVDZ/C\naug-cc-pVDZ/C\ncc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\n"
    reference_basis_block_with_ecp = "aug-cc-pVDZ/C\naug-cc-pVDZ/C\ncc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\naug-cc-pVDZ/C\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\nno-basis-set\n"

    assert generated_basis_block_without_ecp == reference_basis_block_without_ecp
    assert generated_basis_block_with_ecp == reference_basis_block_with_ecp


def test_MRCCInputGenerator_generate_coords_block(mrcc_input_generator):
    mrcc_input_generator_nocp = deepcopy(mrcc_input_generator)

    mrcc_input_generator_nocp.include_cp = False
    mrcc_input_generator_nocp.generate_coords_block()

    mrcc_input_generator.generate_coords_block()

    reference_block_collated = {
        "adsorbate_slab": {
            "float": [
                21.0,
                -2.12018425659,
                -2.12018425659,
                -0.04367284424,
                0.0,
                0.0,
                -0.04269731856,
            ],
            "string": ["charge=-8", "C", "O", "Mg", "Mg", "Mg"],
        },
        "adsorbate": {
            "float": [8.0, -2.12018425659, -2.12018425659],
            "string": ["charge=0", "C", "O"],
        },
        "slab": {
            "float": [
                21.0,
                -2.12018425659,
                -2.12018425659,
                -0.04367284424,
                0.0,
                0.0,
                -0.04269731856,
            ],
            "string": ["charge=-8", "C", "O", "Mg", "Mg", "Mg"],
        },
    }

    reference_block_nocp_collated = {
        "adsorbate_slab": {
            "float": [
                21.0,
                -2.12018425659,
                -2.12018425659,
                -0.04367284424,
                0.0,
                0.0,
                -0.04269731856,
            ],
            "string": ["charge=-8", "C", "O", "Mg", "Mg", "Mg"],
        },
        "adsorbate": {"float": [2.0], "string": ["charge=0", "C"]},
        "slab": {
            "float": [19.0, 2.12018425659, 2.11144262254, -0.04367284424, 0.0, 0.0],
            "string": ["charge=-8", "Mg", "O", "Mg", "Mg"],
        },
    }

    generated_block_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }
    generated_block_nocp_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }

    for system in ["adsorbate_slab", "adsorbate", "slab"]:
        generated_block_collated[system]["float"] = [
            float(x)
            for x in mrcc_input_generator.mrccblocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::10]
        generated_block_collated[system]["string"] = [
            x
            for x in mrcc_input_generator.mrccblocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::5]

        generated_block_nocp_collated[system]["float"] = [
            float(x)
            for x in mrcc_input_generator_nocp.mrccblocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::10]
        generated_block_nocp_collated[system]["string"] = [
            x
            for x in mrcc_input_generator_nocp.mrccblocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::5]

        assert_equal(
            generated_block_collated[system]["string"],
            reference_block_collated[system]["string"],
        )
        assert_allclose(
            generated_block_collated[system]["float"],
            reference_block_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )

        assert_equal(
            generated_block_nocp_collated[system]["string"],
            reference_block_nocp_collated[system]["string"],
        )
        assert_allclose(
            generated_block_nocp_collated[system]["float"],
            reference_block_nocp_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )


def test_MRCCInputGenerator_generate_point_charge_block(mrcc_input_generator):
    generated_point_charge_block = mrcc_input_generator.generate_point_charge_block()

    generated_point_charge_block_shortened = [
        float(x) for x in generated_point_charge_block.split()[5::180]
    ]

    reference_point_charge_block_shortened = [
        -0.04367284424,
        -0.03992370948,
        -2.14923989662,
        -6.37814204923,
        -2.1415520695,
        -4.26789528527,
        -2.1415520695,
        -0.03992370948,
        0.0,
    ]

    assert_allclose(
        generated_point_charge_block_shortened,
        reference_point_charge_block_shortened,
        rtol=1e-05,
        atol=1e-07,
    )


def test_ORCAInputGenerator_init(embedded_adsorbed_cluster, element_info):
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
    orca_input_generator = ORCAInputGenerator(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        ecp_region_indices=[8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info=element_info,
        include_cp=True,
        pal_nprocs_block=pal_nprocs_block,
        method_block=method_block,
        scf_block=scf_block,
        ecp_info=ecp_info,
    )

    # Check when multiplicities is not provided
    assert orca_input_generator.multiplicities == {
        "adsorbate_slab": 1,
        "adsorbate": 1,
        "slab": 1,
    }

    orca_input_generator = ORCAInputGenerator(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
        ecp_region_indices=[8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
        element_info=element_info,
        include_cp=True,
        multiplicities={"adsorbate_slab": 3, "adsorbate": 1, "slab": 2},
        pal_nprocs_block=pal_nprocs_block,
        method_block=method_block,
        scf_block=scf_block,
        ecp_info=ecp_info,
    )

    assert not compare_atoms(
        orca_input_generator.embedded_adsorbed_cluster, embedded_adsorbed_cluster
    )
    assert_equal(orca_input_generator.quantum_cluster_indices, [0, 1, 2, 3, 4, 5, 6, 7])
    assert_equal(orca_input_generator.adsorbate_indices, [0, 1])
    assert_equal(orca_input_generator.slab_indices, [2, 3, 4, 5, 6, 7])
    assert_equal(
        orca_input_generator.ecp_region_indices,
        [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
    )
    assert orca_input_generator.element_info == element_info
    assert orca_input_generator.include_cp is True
    assert orca_input_generator.multiplicities == {
        "adsorbate_slab": 3,
        "adsorbate": 1,
        "slab": 2,
    }

    assert orca_input_generator.pal_nprocs_block == pal_nprocs_block
    assert orca_input_generator.method_block == method_block
    assert orca_input_generator.scf_block == scf_block
    assert orca_input_generator.ecp_info == ecp_info

    # Check if error raise if quantum_cluster_indices and ecp_region_indices overlap

    with pytest.raises(
        ValueError, match="An atom in the quantum cluster is also in the ECP region."
    ):
        orca_input_generator = ORCAInputGenerator(
            embedded_adsorbed_cluster=embedded_adsorbed_cluster,
            quantum_cluster_indices=[0, 1, 2, 3, 4, 5, 6, 7],
            ecp_region_indices=[7, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24],
            element_info=element_info,
            include_cp=True,
            multiplicities={"adsorbate_slab": 3, "adsorbate": 1, "slab": 2},
            pal_nprocs_block=pal_nprocs_block,
            method_block=method_block,
            scf_block=scf_block,
            ecp_info=ecp_info,
        )


def test_ORCAInputGenerator_generate_input(orca_input_generator):
    orca_input_generator_nocp = deepcopy(orca_input_generator)

    orca_input_generator_nocp.include_cp = False
    orca_input_generator_nocp.generate_input()

    orca_input_generator.generate_input()

    reference_block_collated = {
        "adsorbate_slab": {
            "float": [1.0, 2.0, 1.0, 0.0, 2.0],
            "string": [
                "%pal",
                "NewNCore",
                "O",
                '"aug-cc-pVDZ/C"',
                '"orca_svp_start.gbw"',
                "O",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
            ],
        },
        "adsorbate": {
            "float": [1.0],
            "string": [
                "%pal",
                "NewNCore",
                "O",
                '"aug-cc-pVDZ/C"',
                '"orca_svp_start.gbw"',
                "O",
            ],
        },
        "slab": {
            "float": [1.0, 2.0, 1.0, 0.0, 2.0],
            "string": [
                "%pal",
                "NewNCore",
                "O",
                '"aug-cc-pVDZ/C"',
                '"orca_svp_start.gbw"',
                "O:",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
            ],
        },
    }

    reference_block_nocp_collated = {
        "adsorbate_slab": {
            "float": [1.0, 2.0, 1.0, 0.0, 2.0],
            "string": [
                "%pal",
                "NewNCore",
                "O",
                '"aug-cc-pVDZ/C"',
                '"orca_svp_start.gbw"',
                "O",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
            ],
        },
        "adsorbate": {
            "float": [1.0],
            "string": [
                "%pal",
                "NewNCore",
                "O",
                '"aug-cc-pVDZ/C"',
                '"orca_svp_start.gbw"',
                "O",
            ],
        },
        "slab": {
            "float": [1.0, 2.0, 2.10705287155, 0.0, 1.0],
            "string": [
                "%pal",
                "NewNCore",
                "O",
                '"aug-cc-pVDZ/C"',
                '"orca_svp_start.gbw"',
                "O",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "end",
            ],
        },
    }

    generated_block_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }
    generated_block_nocp_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }

    for system in ["adsorbate_slab", "adsorbate", "slab"]:
        generated_block_collated[system]["float"] = [
            float(x)
            for x in orca_input_generator.orcablocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::77]
        generated_block_collated[system]["string"] = [
            x
            for x in orca_input_generator.orcablocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::17]

        assert_equal(
            reference_block_collated[system]["string"],
            generated_block_collated[system]["string"],
        )
        assert_allclose(
            generated_block_collated[system]["float"],
            reference_block_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )

        generated_block_nocp_collated[system]["float"] = [
            float(x)
            for x in orca_input_generator_nocp.orcablocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::77]
        generated_block_nocp_collated[system]["string"] = [
            x
            for x in orca_input_generator_nocp.orcablocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::17]

        assert_equal(
            reference_block_nocp_collated[system]["string"],
            generated_block_nocp_collated[system]["string"],
        )
        assert_allclose(
            generated_block_nocp_collated[system]["float"],
            reference_block_nocp_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )


def test_create_atom_coord_string(embedded_adsorbed_cluster):
    atom = embedded_adsorbed_cluster[0]

    # First let's try the case where it's a normal atom.
    atom_coord_string = create_atom_coord_string(atom=atom)

    with pytest.raises(
        ValueError, match="ECP info cannot be provided for ghost atoms."
    ):
        create_atom_coord_string(
            atom, atom_ecp_info="NewECP\nECP_info1\nECP_info2\n", is_ghost_atom=True
        )

    with pytest.raises(
        ValueError, match="Point charge value must be given for atoms with ECP info."
    ):
        create_atom_coord_string(atom, atom_ecp_info="NewECP\nECP_info1\nECP_info2\n")

    assert (
        atom_coord_string
        == "C                       0.00000000000    0.00000000000    2.00000000000\n"
    )

    # Let's now try the case where it is a ghost atom.
    atom_coord_string = create_atom_coord_string(atom=atom, is_ghost_atom=True)
    assert (
        atom_coord_string
        == "C:                      0.00000000000    0.00000000000    2.00000000000\n"
    )

    # Let's now try the case where it is an atom in the ECP region.
    atom_coord_string = create_atom_coord_string(
        atom=atom, atom_ecp_info="NewECP\nECP_info1\nECP_info2\n", pc_charge=2.0
    )
    assert (
        atom_coord_string
        == "C>     2.00000000000    0.00000000000    0.00000000000    2.00000000000\nNewECP\nECP_info1\nECP_info2\n"
    )


def test_ORCAInputGenerator_generate_coords_block(orca_input_generator):
    orca_input_generator_nocp = deepcopy(orca_input_generator)

    orca_input_generator_nocp.include_cp = False
    orca_input_generator_nocp.generate_coords_block()

    orca_input_generator.generate_coords_block()

    reference_block_collated = {
        "adsorbate_slab": {
            "float": [3.0, 1.0, 5.1757, 1.0, 0.0, 2.0, 1.0],
            "string": [
                "%coords",
                "coords",
                "O",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "end",
            ],
        },
        "adsorbate": {"float": [1.0], "string": ["%coords", "coords", "O:"]},
        "slab": {
            "float": [2.0, 1.0, 5.1757, 1.0, 0.0, 2.0, 1.0],
            "string": [
                "%coords",
                "coords",
                "O",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "end",
            ],
        },
    }

    reference_block_nocp_collated = {
        "adsorbate_slab": {
            "float": [3.0, 1.0, 5.1757, 1.0, 0.0, 2.0, 1.0],
            "string": [
                "%coords",
                "coords",
                "O",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "end",
            ],
        },
        "adsorbate": {"float": [1.0], "string": ["%coords", "coords"]},
        "slab": {
            "float": [2.0, 1.115, 2.0, 2.10705287155, 14.676, 1.0, 1.0],
            "string": [
                "%coords",
                "coords",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
                "p",
                "lmax",
                "Mg>",
                "d",
                "f",
                "NewECP",
                "f",
                "s",
                "N_core",
                "end",
                "p",
            ],
        },
    }

    generated_block_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }
    generated_block_nocp_collated = {
        system: {"float": [], "string": []}
        for system in ["adsorbate_slab", "adsorbate", "slab"]
    }

    for system in ["adsorbate_slab", "adsorbate", "slab"]:
        generated_block_collated[system]["float"] = [
            float(x)
            for x in orca_input_generator.orcablocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::57]
        generated_block_collated[system]["string"] = [
            x
            for x in orca_input_generator.orcablocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::7]

        assert_equal(
            reference_block_collated[system]["string"],
            generated_block_collated[system]["string"],
        )
        assert_allclose(
            generated_block_collated[system]["float"],
            reference_block_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )

        generated_block_nocp_collated[system]["float"] = [
            float(x)
            for x in orca_input_generator_nocp.orcablocks[system].split()
            if x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::57]
        generated_block_nocp_collated[system]["string"] = [
            x
            for x in orca_input_generator_nocp.orcablocks[system].split()
            if not x.replace(".", "", 1).replace("-", "", 1).isdigit()
        ][::7]

        assert_equal(
            reference_block_nocp_collated[system]["string"],
            generated_block_nocp_collated[system]["string"],
        )
        assert_allclose(
            generated_block_nocp_collated[system]["float"],
            reference_block_nocp_collated[system]["float"],
            rtol=1e-05,
            atol=1e-07,
        )


def test_ORCAInputGenerator_format_ecp_info(orca_input_generator):
    with pytest.raises(ValueError):
        orca_input_generator.format_ecp_info(atom_ecp_info="dummy_info\nN_core0\nend")

    atom_ecp_info = """
NewECP
N_core 0
lmax s
s 1
1      1.732000000   14.676000000 2
end
"""
    formatted_atom_ecp_info = orca_input_generator.format_ecp_info(
        atom_ecp_info=atom_ecp_info
    )
    assert (
        formatted_atom_ecp_info
        == "NewECP\nN_core 0\nlmax s\ns 1\n1      1.732000000   14.676000000 2\nend\n"
    )


def test_ORCAInputGenerator_generate_preamble_block(orca_input_generator):
    # Make copy of orca_input_generator for further tests
    orca_input_generator_1 = deepcopy(orca_input_generator)
    orca_input_generator_2 = deepcopy(orca_input_generator)
    orca_input_generator_3 = deepcopy(orca_input_generator)

    # Generate the orca input preamble
    orca_input_generator_1.generate_preamble_block()

    assert (
        orca_input_generator_1.orcablocks["adsorbate_slab"]
        == '%pal nprocs 1 end\n%maxcore 5000 end\n%pointcharges "orca.pc"\n%method\nMethod hf\nRI on\nRunTyp Energy\nNewNCore C 2 end\nNewNCore Mg 2 end\nNewNCore O 2 end\nend\n%basis\nNewGTO C "aug-cc-pVDZ" end\nNewGTO Mg "cc-pVDZ" end\nNewGTO O "aug-cc-pVDZ" end\nNewAuxJGTO C "def2/J" end\nNewAuxJGTO Mg "def2/J" end\nNewAuxJGTO O "def2/JK" end\nNewAuxCGTO C "aug-cc-pVDZ/C" end\nNewAuxCGTO Mg "cc-pVDZ/C" end\nNewAuxCGTO O "aug-cc-pVDZ/C" end\nend\n%scf\nHFTyp rhf\nGuess MORead\nMOInp "orca_svp_start.gbw"\nSCFMode Direct\nsthresh 1e-6\nAutoTRAHIter 60\nMaxIter 1000\nend\n'
    )

    assert (
        orca_input_generator_1.orcablocks["adsorbate_slab"]
        == orca_input_generator_1.orcablocks["adsorbate"]
    )
    assert (
        orca_input_generator_1.orcablocks["adsorbate_slab"]
        == orca_input_generator_1.orcablocks["slab"]
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
    orca_input_generator_2.element_info = element_info
    orca_input_generator_2.generate_preamble_block()

    assert (
        orca_input_generator_2.orcablocks["adsorbate_slab"]
        == '%pal nprocs 1 end\n%maxcore 5000 end\n%pointcharges "orca.pc"\n%method\nMethod hf\nRI on\nRunTyp Energy\nNewNCore C 2 end\nNewNCore Mg 2 end\nNewNCore O 2 end\nend\n%basis\nBasis def2-SVP\nAux def2/J\nAuxC def2-SVP/C\nend\n%scf\nHFTyp rhf\nGuess MORead\nMOInp "orca_svp_start.gbw"\nSCFMode Direct\nsthresh 1e-6\nAutoTRAHIter 60\nMaxIter 1000\nend\n'
    )

    # Testing the case if we provide no blocks
    orca_input_generator_3.scf_block = None
    orca_input_generator_3.method_block = None
    orca_input_generator_3.pal_nprocs_block = None
    orca_input_generator_3.element_info = None
    orca_input_generator_3.generate_preamble_block()

    assert (
        orca_input_generator_3.orcablocks["adsorbate_slab"]
        == '%pointcharges "orca.pc"\n'
    )

    # Check whether error raised if not all element_info is provided
    with pytest.raises(ValueError):
        element_info_error = {"C": element_info["C"]}
        orca_input_generator_3.element_info = element_info_error
        orca_input_generator_3.generate_preamble_block()


def test_ORCAInputGenerator_create_point_charge_file(orca_input_generator, tmp_path):
    # Create the point charge file
    orca_input_generator.create_point_charge_file(pc_file=tmp_path / "orca.pc")

    # Read the written file
    orca_pc_file = np.loadtxt(tmp_path / "orca.pc", skiprows=1)

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
        adsorbate_slab_file=Path(FILE_DIR, "skzcam_files", "NO_MgO.poscar.gz"),
        slab_center_indices=[32, 33],
        adsorbate_indices=[0, 1],
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
    assert_equal(adsorbate.get_atomic_numbers().tolist(), [7, 8])

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
    assert_equal(
        slab.get_atomic_numbers().tolist(),
        [
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
        ],
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
        slab=slab,
        slab_center_idx=30,
        atom_oxi_states={"Mg": 2.0, "O": -2.0},
        filepath=tmp_path,
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

    assert_equal(
        chemshell_embedded_cluster.get_atomic_numbers()[::20].tolist(),
        [
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
        ],
    )


def test_convert_pun_to_atoms():
    embedded_cluster = convert_pun_to_atoms(
        pun_file=Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
        atom_oxi_states={"Mg": 2.0, "O": -2.0},
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
    assert_equal(
        embedded_cluster.get_array("atom_type")[:10],
        [
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
        ],
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
        embedded_cluster=embedded_cluster,
        adsorbate=adsorbate,
        adsorbate_vector_from_slab=[0.0, 0.0, 2.0],
        quantum_cluster_indices=[[0, 1, 3, 4], [5, 6, 7, 8]],
        ecp_region_indices=[[0, 1, 3, 4], [5, 6, 7, 8]],
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

    assert_equal(
        embedded_cluster.get_chemical_symbols()[:10],
        ["C", "O", "Mg", "O", "O", "O", "O", "O", "Mg", "Mg"],
    )
    assert_allclose(
        embedded_cluster.get_array("oxi_states")[:10],
        np.array([0.0, 0.0, 2.0, -2.0, -2.0, -2.0, -2.0, -2.0, 2.0, 2.0]),
        rtol=1e-05,
        atol=1e-07,
    )
    assert_equal(
        embedded_cluster.get_array("atom_type")[:10],
        [
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
        ],
    )

    # Check that the quantum_idx and ecp_idx match the reference
    assert_equal(quantum_idx, [[0, 1, 2, 3, 5, 6], [0, 1, 7, 8, 9, 10]])
    assert_equal(ecp_idx, [[2, 3, 5, 6], [7, 8, 9, 10]])


def test_get_atom_distances():
    # Creating a H2 molecule as an Atoms object
    h2_molecule = Atoms("H2", positions=[(0, 0, 0), (0, 0, 2)])

    # Run _get_atom_distances function to get distance of h2 molecule atoms from a center position
    atom_distances = _get_atom_distances(
        embedded_cluster=h2_molecule, center_position=[2, 0, 0]
    )

    assert_allclose(atom_distances, np.array([2.0, 2.82842712]), rtol=1e-05, atol=1e-07)


def test_find_cation_shells(embedded_cluster):
    # Get distance of atoms from the center
    distances = _get_atom_distances(
        embedded_cluster=embedded_cluster, center_position=[0, 0, 2]
    )

    # Find the cation shells from the distances
    cation_shells, cation_shells_idx = _find_cation_shells(
        embedded_cluster=embedded_cluster, distances=distances, shell_width=0.005
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

    assert_equal(
        cation_shells_idx_flatten, [0, 9, 8, 6, 7, 11, 12, 10, 13, 19, 21, 18, 20, 22]
    )


def test_get_anion_coordination(embedded_cluster, distance_matrix):
    # Get the anions for the second SKZCAM shell
    anion_shell_idx = _get_anion_coordination(
        embedded_cluster=embedded_cluster,
        cation_shell_indices=[9, 8, 6, 7],
        dist_matrix=distance_matrix,
    )

    # Check anion indices match with reference
    assert_equal(
        anion_shell_idx, [1, 2, 3, 4, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 29, 30]
    )


def test_get_ecp_region(embedded_cluster, distance_matrix):
    # Find the ECP region for the first cluster
    ecp_region_idx = _get_ecp_region(
        embedded_cluster=embedded_cluster,
        quantum_cluster_indices=[[0, 1, 2, 3, 4, 5]],
        dist_matrix=distance_matrix,
        ecp_dist=3,
    )

    # Check ECP region indices match with reference
    assert_equal(ecp_region_idx[0], [6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22])


def test_create_skzcam_clusters(tmp_path):
    # Get quantum cluster and ECP region indices
    _, quantum_cluster_idx, ecp_region_idx = create_skzcam_clusters(
        pun_file=Path(FILE_DIR, "skzcam_files", "mgo_shells_cluster.pun.gz"),
        center_position=[0, 0, 2],
        atom_oxi_states={"Mg": 2.0, "O": -2.0},
        shell_max=2,
        ecp_dist=3.0,
        shell_width=0.005,
        write_clusters=True,
        write_clusters_path=tmp_path,
    )

    # Check quantum cluster indices match with reference
    assert_equal(
        quantum_cluster_idx[1],
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 29, 30],
    )

    # Check ECP region indices match with reference
    assert_equal(
        ecp_region_idx[1],
        [
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
        ],
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

    assert_equal(skzcam_cluster.get_atomic_numbers().tolist(), [12, 8, 8, 8, 8, 8])
