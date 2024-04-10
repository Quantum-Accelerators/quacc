from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, fcc100, molecule
from ase.io import read
from ase import Atoms

from quacc.atoms.skzcam import (
	create_skzcam_clusters,
	convert_pun_to_atoms,
	_get_atom_distances,
	_find_cation_shells,
	_get_anion_coordination,
	_get_ecp_region
)

FILE_DIR = Path(__file__).parent

@pytest.fixture
def embedded_cluster():
    return convert_pun_to_atoms( FILE_DIR / "mgo_shells_cluster.pun.gz", {'Mg': 2.0, 'O': -2.0})

@pytest.fixture
def distance_matrix(embedded_cluster):
    return embedded_cluster.get_all_distances()

def test_convert_pun_to_atoms():
	embedded_cluster = convert_pun_to_atoms( FILE_DIR / "mgo_shells_cluster.pun.gz", {'Mg': 2.0, 'O': -2.0})
	
	# Check that number of atoms matches our reference
	assert len(embedded_cluster) == 390

	# Check that last 10 elements of the oxi_state match our reference
	np.testing.assert_allclose(embedded_cluster.get_array("oxi_states")[-10:], np.array([-0.80812511, 2.14427889, -0.96000248, 2.14427887, -0.8081251, 2.10472993, -0.89052904,2.10472993, -0.8081251, 2.14427887]), rtol=1e-05, atol=1e-07)
	
	# Check that first 10 elements of atom_type array match our reference
	assert np.all(embedded_cluster.get_array("atom_type")[:10] == ['cation', 'anion', 'anion', 'anion' ,'anion', 'anion', 'cation', 'cation', 'cation', 'cation'])
	
	# Check that the positions of the atom matches
	np.testing.assert_allclose(embedded_cluster[200].position, np.array([ 6.33074029, -2.11024676, -6.37814205]), rtol=1e-05, atol=1e-07) 
	
def test_get_atom_distances():
	
	# Creating a H2 molecule as an Atoms object
	h2_molecule = Atoms('H2', positions=[(0, 0, 0), (0, 0, 2)])
	
	# Run _get_atom_distances function to get distance of h2 molecule atoms from a center position
	atom_distances = _get_atom_distances(h2_molecule,[2,0,0])
	
	np.testing.assert_allclose(atom_distances,np.array([2.        , 2.82842712]),rtol=1e-05, atol=1e-07)	

def test_find_cation_shells(embedded_cluster):
	
	# Get distance of atoms from the center
	distances = _get_atom_distances(embedded_cluster,[0,0,2])

	# Find the cation shells from the distances
	cation_shells, cation_shells_idx = _find_cation_shells(embedded_cluster,distances)

	# As these list of lists do not have the same length, we flatten first 5 lists into a 1D list for comparison
	cation_shells_flatten = [item for row in cation_shells[:5] for item in row] 
	cation_shells_idx_flatten = [item for row in cation_shells_idx[:5] for item in row]

	# Check that these lists are correct
	np.testing.assert_allclose(cation_shells_flatten,np.array([2.0,
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
 6.267895285274443]),rtol=1e-05, atol=1e-07)
	
	assert np.all(cation_shells_idx_flatten == [0, 9, 8, 6, 7, 11, 12, 10, 13, 19, 21, 18, 20, 22])

def test_get_anion_coordination(embedded_cluster, distance_matrix):

	# Get the anions for the second SKZCAM shell
	anion_shell_idx = _get_anion_coordination(embedded_cluster, [9, 8, 6, 7], distance_matrix)

	# Check anion indices match with reference
	assert np.all(anion_shell_idx == [1, 2, 3, 4, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 29, 30])


def test_get_ecp_region(embedded_cluster, distance_matrix):
	
	# Find the ECP region for the first cluster
	ecp_region_idx = _get_ecp_region(embedded_cluster, quantum_cluster_idx = [[0,1,2,3,4,5]], dist_matrix = distance_matrix, ecp_dist = 3)

	# Check ECP region indices match with reference
	assert np.all(ecp_region_idx[0] == [6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22])
	
def test_create_skzcam_clusters():

	# Get quantum cluster and ECP region indices
	quantum_cluster_idx, ecp_region_idx = create_skzcam_clusters(FILE_DIR / "mgo_shells_cluster.pun.gz", [0,0,2], {'Mg': 2.0, 'O': -2.0}, shell_max=2,ecp_dist=3.0)

	# Check quantum cluster indices match with reference
	assert np.all(quantum_cluster_idx[1] == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 29, 30])


	# Check ECP region indices match with reference
	assert np.all(ecp_region_idx[1] == [10,
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
 83]
)
