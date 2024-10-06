from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_numbers
from ase.io import read, write
from ase.io.orca import write_orca
from ase.units import Bohr
from monty.dev import requires
from monty.io import zopen
from monty.os.path import zpath

from quacc.calculators.mrcc.io import write_mrcc

if TYPE_CHECKING:
    from ase.atom import Atom
    from numpy.typing import NDArray

    from quacc.types import (
        BlockInfo,
        ElementInfo,
        ElementStr,
        MRCCInputDict,
        MultiplicityDict,
        SKZCAMInfo,
        SKZCAMOutput,
        semicore_frozen_core_dict,
        valence_frozen_core_dict,
    )


element_from_atomic_num_dict = {v: k for k, v in atomic_numbers.items()}

semicore_frozen_core_dict = {
    element_from_atomic_num_dict[x]: (
        0 if x <= 4 else 
        2 if x <= 12 else 
        10 if x <= 30 else 
        18 if x <= 38 else 
        28 if x <= 48 else 
        36 if x <= 71 else 
        46 if x <= 80 else 
        68 if x <= 103 else 
        None  # You can choose a default value or None for x > 103
    )
    for x in range(1, 104)  # Adjust the range to include up to 103
    if x in element_from_atomic_num_dict  # Ensure the atomic number exists in the dictionary
}

valence_frozen_core_dict = {
    element_from_atomic_num_dict[x]: (
        0 if x <= 2 else 
        2 if x <= 10 else 
        10 if x <= 18 else 
        18 if x <= 30 else 
        28 if x <= 36 else 
        36 if x <= 48 else 
        46 if x <= 54 else 
        54 if x <= 70 else 
        68 if x <= 80 else 
        78 if x <= 86 else 
        86 if x <= 102 else
        100 if x <= 103 else
        None  # You can choose a default value or None for x > 103
    )
    for x in range(1, 104)  # Adjust the range to include up to 103
    if x in element_from_atomic_num_dict  # Ensure the atomic number exists in the dictionary
}

skzcam_cation_cap_ecp = {
    'Ti': """NewECP
N_core 0
  lmax f
  s 2
   1      0.860000       9.191690  2
   2      0.172000       0.008301  2
  p 2
   1      0.860000      17.997720  2
   2      0.172000      -0.032600  2
  d 2
   1      1.600000      -9.504310  2
   2      0.320000      -0.151370  2
  f 1
   1      1.000000000    0.000000000 2
end""",
    'Mg': """NewECP
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

has_chemshell = find_spec("chemsh") is not None

class SKZCAMInputSet:
    """
    A class to generate the complete set of SKZCAM inputs for performing direct MRCC and/or ORCA calculations. We define different ONIOM """

    def __init__(
            self,
            adsorbate_slab_embedded_cluster: Atoms,
            quantum_cluster_indices_set: list[list[int]],
            ecp_region_indices_set: list[list[int]],
            mp2_oniom1_ll: SKZCAMInfo | None = None,
            mp2_oniom1_hl: SKZCAMInfo | None = None,
            mp2_oniom2_hl: SKZCAMInfo | None = None,
            mp2_oniom3_hl: SKZCAMInfo | None = None,
            mp2_oniom4_hl: SKZCAMInfo | None = None,
            deltaCC: SKZCAMInfo | None = None 
    ):
        """
        Parameters
        ----------
        adsorbate_slab_embedded_cluster
            The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        quantum_cluster_indices_set
            A list of lists containing the indices of the atoms of a set of quantum clusters. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        ecp_region_indices_set
            A list of lists containing the indices of the atoms in the ECP region of a set of quantum clusters. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        mp2_oniom1_ll
            A dictionary containing the information for the "low-level" MP2 treatment of ONIOM layer 1.
        mp2_oniom1_hl
            A dictionary containing the information for the "high-level" MP2 treatment of ONIOM layer 1.
        mp2_oniom2_hl
            A dictionary containing the information for the "high-level" MP2 treatment of ONIOM layer 2.
        mp2_oniom3_hl
            A dictionary containing the information for the "high-level" MP2 treatment of ONIOM layer 3.
        mp2_oniom4_hl
            A dictionary containing the information for the "high-level" MP2 treatment of ONIOM layer 4.
        deltaCC
            A dictionary containing the information for the deltaCC calculation.

        Returns
        -------
        None
        """

        def _format_oniom_dict(
                oniom_dict: SKZCAMInfo
        ) -> SKZCAMInfo:
            """
            Format the oniom dictionary to include the default values for the SKZCAM calculations if the user hasn't specified it themselves.

            Parameters
            ----------
            oniom_dict
                The dictionary containing the information for the SKZCAM calculation.

            Returns
            -------
            SKZCAMInfo
                The formatted dictionary containing the information for the SKZCAM calculation.
            
            """

            # Raise error if the code is not MRCC or ORCA
            if oniom_dict['code'] not in ['mrcc', 'orca']:
                raise ValueError("The code must be either 'mrcc' or 'orca'.")
            
            if 'multiplicities' not in oniom_dict or oniom_dict['multiplicities'] is None:
                oniom_dict['multiplicities'] = {"adsorbate_slab": 1, "adsorbate": 1, "slab": 1}
            if 'ecp' not in oniom_dict or oniom_dict['ecp'] is None:
                oniom_dict['ecp'] = {}
            if 'ri_scf_basis' not in oniom_dict:
                oniom_dict['ri_scf_basis'] = None
            if 'ri_cwft_basis' not in oniom_dict:
                oniom_dict['ri_cwft_basis'] = None
            if 'nprocs' not in oniom_dict or oniom_dict['nprocs'] is None:
                oniom_dict['nprocs'] = 1
            if 'max_memory' not in oniom_dict or oniom_dict['max_memory'] is None:
                oniom_dict['max_memory'] = 1000
            if 'mrcc_calc_inputs' not in oniom_dict or oniom_dict['mrcc_calc_inputs'] is None:
                oniom_dict['mrcc_calc_inputs'] = {}

            if oniom_dict['code'] == 'orca':
                if 'orca_method_block' not in oniom_dict or oniom_dict['orca_method_block'] is None:
                    oniom_dict['orca_method_block'] = {"RI": "on",
                                                       "RunTyp": "Energy"}
                if 'orca_scf_block' not in oniom_dict or oniom_dict['orca_scf_block'] is None:
                    oniom_dict['orca_scf_block'] = {"HFTyp": "rhf",
                                                    "Guess": "PAtom",
                                                    "SCFMode": "Direct",
                                                    "sthresh": "1e-6",
                                                    "AutoTRAHIter": "60",
                                                    "MaxIter": "1000"}
            if 'orca_cation_cap_ecp' not in oniom_dict or oniom_dict['orca_cation_cap_ecp'] is None:
                oniom_dict['orca_cation_cap_ecp'] = skzcam_cation_cap_ecp
            return oniom_dict
                
        self.adsorbate_slab_embedded_cluster = adsorbate_slab_embedded_cluster
        self.quantum_cluster_indices_set = quantum_cluster_indices_set
        self.ecp_region_indices_set = ecp_region_indices_set
        self.skzcam_input_sets = {}
        if mp2_oniom1_ll:
            self.mp2_oniom1_ll = mp2_oniom1_ll
            self.skzcam_input_sets["mp2_oniom1_ll"] = _format_oniom_dict(mp2_oniom1_ll)
        if mp2_oniom1_hl:
            self.mp2_oniom1_hl = mp2_oniom1_hl
            self.skzcam_input_sets["mp2_oniom1_hl"] = _format_oniom_dict(mp2_oniom1_hl)
        if mp2_oniom2_hl:
            self.mp2_oniom2_hl = mp2_oniom2_hl
            self.skzcam_input_sets["mp2_oniom2_hl"] = _format_oniom_dict(mp2_oniom2_hl)
        if mp2_oniom3_hl:
            self.mp2_oniom3_hl = mp2_oniom3_hl
            self.skzcam_input_sets["mp2_oniom3_hl"] = _format_oniom_dict(mp2_oniom3_hl)
        if mp2_oniom4_hl:
            self.mp2_oniom4_hl = mp2_oniom4_hl
            self.skzcam_input_sets["mp2_oniom4_hl"] = _format_oniom_dict(mp2_oniom4_hl)
        if deltaCC:
            self.deltaCC = deltaCC
            self.skzcam_input_sets["deltaCC"] = deltaCC


    def create_element_info(
            self,
            frozencore: Literal['valence', 'semicore'] | dict[ElementStr, int],
            basis : Literal['DZ', 'TZ', 'QZ', '5Z'] | dict[ElementStr, str],
            code: Literal['mrcc', 'orca'],
            ecp: dict[ElementStr, str] | None = None,
            ri_scf_basis: dict[ElementStr, str] | None = None,
            ri_cwft_basis: dict[ElementStr, str] | None = None
    ):
        """
        Creates the element info dictionary for the SKZCAM input across each oniom layer.

        Parameters
        ----------
        basis
            The basis set to use for the quantum cluster. This could be either double-zeta, triple-zeta, quadruple-zeta, quintuple-zeta, denoted as 'DZ', 'TZ',' QZ' and '5Z' respectively.
        frozencore
            The frozen core to use for the quantum cluster. This could be specified as a string being either 'semicore' or 'valence'
        code
            The code to use for the quantum cluster. This could be either mrcc or orca.
        ecp
            The effective core potential to use for each element within the quantum cluster.

        Returns
        -------
        dict[ElementStr, ElementInfo]
            A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
        """

        # Create an adsorbate_slab_quantum_cluster object for the first SKZCAM cluster
        adsorbate_slab_quantum_cluster = self.adsorbate_slab_embedded_cluster[self.quantum_cluster_indices_set[0]]
        
        # Define set of elements
        set_of_elements = list(set(adsorbate_slab_quantum_cluster.get_chemical_symbols()))

        # Choose whether to provide preset values for the element info dictionary based on the SKZCAM paper:
        # Preset values are used if basis is not a dictionary or frozencore is not a dictionary
        use_presets = False
        if isinstance(basis, str) and isinstance(frozencore, str):
            if basis in ['DZ','TZ','QZ','5Z'] and frozencore in ['valence','semicore']:
                use_presets = True

        element_info_dict = {}
        # If use_presets is True, use some preset inputs based on basis set and frozen core
        if use_presets:
            for atom_idx, atom in enumerate(adsorbate_slab_quantum_cluster):
                if atom.symbol in element_info_dict:
                    continue
                if adsorbate_slab_quantum_cluster.get_array("atom_type")[atom_idx] in ['adsorbate','anion']:
                    element_info_dict[atom.symbol] =  {
                        'core' : valence_frozen_core_dict[atom.symbol],
                        'basis' : f'aug-cc-pV{basis}',
                        'ecp': ecp[atom.symbol] if atom.symbol in ecp else 'none',
                        'ri_scf_basis': 'def2-QZVPP-RI-JK' if code == 'mrcc' else 'def2/J',
                        'ri_cwft_basis': f'aug-cc-pV{basis}-RI' if code == 'mrcc' else f'aug-cc-pV{basis}/C'
                    }
                elif adsorbate_slab_quantum_cluster.get_array("atom_type")[atom_idx] == 'cation' and frozencore == 'valence':
                    element_info_dict[atom.symbol] =  {
                        'core' : valence_frozen_core_dict[atom.symbol],
                        'basis' : f'cc-pV{basis}',
                        'ecp': ecp[atom.symbol] if atom.symbol in ecp else 'none',
                        'ri_scf_basis': 'def2-QZVPP-RI-JK' if code == 'mrcc' else 'def2/J',
                        'ri_cwft_basis': f'cc-pV{basis}-RI' if code == 'mrcc' else f'cc-pV{basis}/C'
                    }
                elif adsorbate_slab_quantum_cluster.get_array("atom_type")[atom_idx] == 'cation' and frozencore == 'semicore':
                    element_info_dict[atom.symbol] =  {
                        'core' : semicore_frozen_core_dict[atom.symbol],
                        'basis' : f'cc-pwCV{basis}',
                        'ecp': ecp[atom.symbol] if atom.symbol in ecp else 'none',
                        'ri_scf_basis': 'def2-QZVPP-RI-JK' if code == 'mrcc' else 'def2/J',
                        'ri_cwft_basis': f'AutoAux' if code == 'mrcc' else f'cc-pwCV{basis}/C'
                    }
        # Use custom values for the element info dictionary
        else:
            # First generate list of elements
            for atom_idx, atom in enumerate(adsorbate_slab_quantum_cluster):
                if atom.symbol in element_info_dict:
                    continue
                element_info_dict[atom.symbol] =  {}
                # Update the element info dictionary with the custom values

                # Firstly for frozen core
                if isinstance(frozencore, dict):
                    # Raise error if frozencore is not provided for all elements
                    if not all([element in frozencore for element in set_of_elements]):
                        raise ValueError("frozencore must be provided for all elements in the quantum cluster when provided as a dictionary.")
                    element_info_dict[atom.symbol]['core'] = frozencore[atom.symbol]
                elif isinstance(frozencore, str):
                    if frozencore not in ['valence','semicore']:
                        raise ValueError("frozencore must be provided as either 'valence' or 'semicore' if provided as a string.")
                    if frozencore == 'valence':
                        element_info_dict[atom.symbol]['core'] = valence_frozen_core_dict[atom.symbol]
                    elif frozencore == 'semicore':
                        element_info_dict[atom.symbol]['core'] = semicore_frozen_core_dict[atom.symbol]
                else:
                    raise ValueError("frozencore must be provided as a string or as a dictionary.")
                
                # Next for basis
                if isinstance(basis, dict):
                    # Raise error if basis is not provided for all elements
                    if not all([element in basis for element in set_of_elements]):
                        raise ValueError("basis must be provided for all elements in the quantum cluster when provided as a dictionary.")
                    element_info_dict[atom.symbol]['basis'] = basis[atom.symbol]
                elif isinstance(basis, str):
                    element_info_dict[atom.symbol]['basis'] = basis
                else:
                    raise ValueError("basis must be provided as a string or dictionary of elements.")
                
                # Next for ecp
                if atom.symbol in ecp:
                    element_info_dict[atom.symbol]['ecp'] = ecp[atom.symbol]
                else:
                    element_info_dict[atom.symbol]['ecp'] = 'none'


                # Next for ri_scf_basis
                if ri_scf_basis is None:
                    element_info_dict[atom.symbol]['ri_scf_basis'] = 'def2-QZVPP-RI-JK' if code == 'mrcc' else 'def2/J'
                elif isinstance(ri_scf_basis, dict):
                    # Raise error if ri_scf_basis is not provided for all elements
                    if not all([element in ri_scf_basis for element in set_of_elements]):
                        raise ValueError("ri_scf_basis must be provided for all elements in the quantum cluster when provided as a dictionary.")
                    element_info_dict[atom.symbol]['ri_scf_basis'] = ri_scf_basis[atom.symbol]
                elif isinstance(ri_scf_basis, str):
                    element_info_dict[atom.symbol]['ri_scf_basis'] = ri_scf_basis
                else:
                    raise ValueError("ri_scf_basis must be provided as a string or dictionary of elements.")


                # Next for ri_cwft_basis
                if ri_cwft_basis is None:
                    element_info_dict[atom.symbol]['ri_cwft_basis'] = element_info_dict[atom.symbol]['basis'] + '-RI' if code == 'mrcc' else element_info_dict[atom.symbol]['basis'] + '/C' 
                elif isinstance(ri_cwft_basis, dict):
                    # Raise error if ri_cwft_basis is not provided for all elements
                    if not all([element in ri_cwft_basis for element in set_of_elements]):
                        raise ValueError("ri_cwft_basis must be provided for all elements in the quantum cluster when provided as a dictionary.")
                    element_info_dict[atom.symbol]['ri_cwft_basis'] = ri_cwft_basis[atom.symbol]
                elif isinstance(ri_cwft_basis, str):
                    element_info_dict[atom.symbol]['ri_cwft_basis'] = ri_cwft_basis
                else:
                    raise ValueError("ri_cwft_basis must be provided as a string or dictionary of elements.")

        return element_info_dict
    
    def generate_input(
            self,
            input_dir: str | Path,
    ):
        """
        Generates the SKZCAM input for the MRCC and ORCA ASE calculators.
        
        Parameters
        ----------
        input_dir
            The directory where the input files will be written.
            
        Returns
        -------
        None
        """

        # Start by writing the input files for the MP2 calculations
        for oniom_method, oniom_parameters in self.skzcam_input_sets.items():
            for cluster_num in range(1,oniom_parameters["max_cluster_num"]+1):
                if isinstance(oniom_parameters['basis'], dict):
                    # Get a list of basis sets for the oniom_method. This would just be the single basis set specified if there is no CBS(../..) specified, otherwise, it would be a list of two basis sets.
                    if all(['CBS(/)' == basis.split('/')[0][:4] + '/' + basis.split('/')[-1][-1:] for basis in oniom_parameters['basis'].values()]):                     
                        basis_sets = [{element: basis[4:-1].split('/')[0] for element, basis in oniom_parameters['basis'].items()}, {element: basis[4:-1].split('/')[1] for element, basis in oniom_parameters['basis'].items()}]
                    else:
                        basis_sets = [oniom_parameters["basis"]]
                elif isinstance(oniom_parameters['basis'], str):
                    if 'CBS(/)' == oniom_parameters['basis'].split('/')[0][:4] + '/' + oniom_parameters['basis'].split('/')[-1][-1:] :
                        basis_sets = [oniom_parameters['basis'][4:-1].split('/')[0], oniom_parameters['basis'][4:-1].split('/')[1]]
                    else:
                        basis_sets = [oniom_parameters["basis"]]
                else:
                    raise ValueError("basis must be given as a string or a dictionary.")
                
                for basis_set in basis_sets:
                    if isinstance(basis_set, dict):
                        ordered_basis_set = dict(sorted(basis_set.items()))
                        basis_name = '_'.join([key + '-' +  value for key, value in ordered_basis_set.items()])
                    elif isinstance(basis_set, str):
                        if isinstance(oniom_parameters['frozencore'], str):
                            if oniom_parameters['frozencore'] == 'valence':
                                basis_name = 'aV' + basis_set
                            elif oniom_parameters['frozencore'] == 'semicore':
                                basis_name = 'awCV' + basis_set
                        else:
                            basis_name = basis_set
                    if oniom_parameters["code"] == 'mrcc':
                        element_info = self.create_element_info(
                            basis=basis_set,
                            frozencore=oniom_parameters["frozencore"],
                            code=oniom_parameters["code"],
                            ecp=oniom_parameters["ecp"],
                            ri_scf_basis=oniom_parameters["ri_scf_basis"],
                            ri_cwft_basis=oniom_parameters["ri_cwft_basis"]
                        )
                        mrcc_block_inputs = MRCCInputGenerator(
                                    adsorbate_slab_embedded_cluster=self.adsorbate_slab_embedded_cluster,
                                    quantum_cluster_indices=self.quantum_cluster_indices_set[cluster_num-1],
                                    ecp_region_indices=self.ecp_region_indices_set[cluster_num-1],
                                    element_info=element_info,
                                    include_cp=True,
                                    multiplicities={"adsorbate_slab": 1, "adsorbate": 1, "slab": 1},
                                ).generate_input()
                        
                        mrcc_default_calc_inputs = {'calc': 'DF-MP2', 'scftype': 'rhf', 'verbosity': 3, 'mem': f"{oniom_parameters['max_memory']}MB", 'symm': 'off', 'unit': 'angs', 'scfiguess': 'small', 'scfmaxit': 1000, 'scfalg': 'locfit1'}

                        mrcc_calc_inputs = {**mrcc_default_calc_inputs,**oniom_parameters['mrcc_calc_inputs']}

                        mrcc_inputs = {structure: {**mrcc_calc_inputs, **mrcc_block_inputs[structure]}  for structure in mrcc_block_inputs}

                        write_mrcc(Path(input_dir,f'MRCC_MINP_MP2_cluster_{cluster_num}_{basis_name}_adsorbate_slab'),self.adsorbate_slab_embedded_cluster,mrcc_inputs['adsorbate_slab'])
                        write_mrcc(Path(input_dir,f'MRCC_MINP_MP2_cluster_{cluster_num}_{basis_name}_slab'),self.adsorbate_slab_embedded_cluster,mrcc_inputs['slab'])
                        write_mrcc(Path(input_dir,f'MRCC_MINP_MP2_cluster_{cluster_num}_{basis_name}_adsorbate'),self.adsorbate_slab_embedded_cluster,mrcc_inputs['adsorbate'])
                    elif oniom_parameters["code"] == 'orca':
                        element_info = self.create_element_info(
                            basis=oniom_parameters["basis"],
                            frozencore=oniom_parameters["frozencore"],
                            code=oniom_parameters["code"],
                            ecp=oniom_parameters["ecp"],
                            ri_scf_basis=oniom_parameters["ri_scf_basis"],
                            ri_cwft_basis=oniom_parameters["ri_cwft_basis"]
                        )

                        orca_input_generator = ORCAInputGenerator(
                            adsorbate_slab_embedded_cluster=self.adsorbate_slab_embedded_cluster,
                            quantum_cluster_indices=self.quantum_cluster_indices_set[cluster_num-1],
                            ecp_region_indices=self.ecp_region_indices_set[cluster_num-1],
                            element_info=element_info,
                            include_cp=True,
                            multiplicities={"adsorbate_slab": 1, "adsorbate": 1, "slab": 1},
                            pal_nprocs_block={"nprocs": int(oniom_parameters["nprocs"]), "maxcore": int(float(oniom_parameters["max_memory"])/int(oniom_parameters["nprocs"]))},
                            method_block=oniom_parameters["orca_method_block"],
                            scf_block=oniom_parameters["orca_scf_block"],
                            ecp_info=oniom_parameters["orca_cation_cap_ecp"]
                        )
                        orca_blocks = orca_input_generator.generate_input()

                        orca_inputs = {structure: {'orcasimpleinput': 'TightSCF RI-MP2 RIJCOSX SlowConv DIIS', 'orcablocks': orca_blocks[structure]}for structure in orca_blocks}

                        write_orca(Path(input_dir,f'ORCA_MP2_cluster_{cluster_num}_{basis_name}_adsorbate.inp'),self.adsorbate_slab_embedded_cluster,orca_inputs['adsorbate'])
                        write_orca(Path(input_dir,f'ORCA_MP2_cluster_{cluster_num}_{basis_name}_slab.inp'),self.adsorbate_slab_embedded_cluster,orca_inputs['slab'])
                        write_orca(Path(input_dir,f'ORCA_MP2_cluster_{cluster_num}_{basis_name}_adsorbate_slab.inp'),self.adsorbate_slab_embedded_cluster,orca_inputs['adsorbate_slab'])
                        orca_input_generator.create_point_charge_file(Path(input_dir,f'ORCA_MP2_cluster_{cluster_num}_{basis_name}.pc'))

            


class MRCCInputGenerator:
    """
    A class to generate the SKZCAM input for the MRCC ASE calculator.

    Attributes
    ----------
    adsorbate_slab_embedded_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created within the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
    quantum_cluster_indices
        A list containing the indices of the atoms in one quantum cluster. These indices are created within the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
    ecp_region_indices
        A list containing the indices of the atoms in one ECP region. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    include_cp
        If True, the coords strings will include the counterpoise correction (i.e., ghost atoms) for the adsorbate and slab.
    multiplicities
        The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.
    adsorbate_slab_cluster
        The ASE Atoms object for the quantum cluster of the adsorbate-slab complex.
    ecp_region
        The ASE Atoms object for the ECP region.
    adsorbate_indices
        The indices of the adsorbates from the adsorbate_slab_cluster quantum cluster.
    slab_indices
        The indices of the slab from the adsorbate_slab_cluster quantum cluster.
        The ECP region cluster.
    adsorbate_cluster
        The ASE Atoms object for the quantum cluster of the adsorbate.
    slab_cluster
        The ASE Atoms object for the quantum cluster of the slab.
    skzcam_input_str
        The MRCC input block (to be put in 'skzcam_input_str' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
    """

    def __init__(
        self,
        adsorbate_slab_embedded_cluster: Atoms,
        quantum_cluster_indices: list[int],
        ecp_region_indices: list[int],
        element_info: dict[ElementStr, ElementInfo] | None = None,
        include_cp: bool = True,
        multiplicities: MultiplicityDict | None = None,
    ) -> None:
        """
        Parameters
        ----------
        adsorbate_slab_embedded_cluster
            The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created within the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        quantum_cluster_indices
            A list containing the indices of the atoms in one quantum cluster. These indices are created within the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        ecp_region_indices
            A list containing the indices of the atoms in the corresponding ECP region of one quantum cluster. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        element_info
            A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
        include_cp
            If True, the coords strings will include the counterpoise correction (i.e., ghost atoms) for the adsorbate and slab.
        multiplicities
            The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.

        Returns
        -------
        None
        """

        self.adsorbate_slab_embedded_cluster = adsorbate_slab_embedded_cluster
        self.quantum_cluster_indices = quantum_cluster_indices
        self.ecp_region_indices = ecp_region_indices
        self.element_info = element_info
        self.include_cp = include_cp
        self.multiplicities = (
            {"adsorbate_slab": 1, "adsorbate": 1, "slab": 1}
            if multiplicities is None
            else multiplicities
        )

        # Check that none of the indices in quantum_cluster_indices are in ecp_region_indices
        if not np.all(
            [x not in self.ecp_region_indices for x in self.quantum_cluster_indices]
        ):
            raise ValueError(
                "An atom in the quantum cluster is also in the ECP region."
            )

        # Create the adsorbate-slab complex quantum cluster and ECP region cluster
        self.adsorbate_slab_cluster: Atoms = self.adsorbate_slab_embedded_cluster[
            self.quantum_cluster_indices
        ]
        self.ecp_region: Atoms = self.adsorbate_slab_embedded_cluster[
            self.ecp_region_indices
        ]

        # Get the indices of the adsorbates from the quantum cluster
        self.adsorbate_indices: list[int] = [
            i
            for i in range(len(self.adsorbate_slab_cluster))
            if self.adsorbate_slab_cluster.get_array("atom_type")[i] == "adsorbate"
        ]
        # Get the indices of the slab from the quantum cluster
        self.slab_indices: list[int] = [
            i
            for i in range(len(self.adsorbate_slab_cluster))
            if self.adsorbate_slab_cluster.get_array("atom_type")[i] != "adsorbate"
        ]

        # Create the adsorbate and slab quantum clusters
        self.adsorbate_cluster: Atoms = self.adsorbate_slab_cluster[
            self.adsorbate_indices
        ]
        self.slab_cluster: Atoms = self.adsorbate_slab_cluster[self.slab_indices]

        # Initialize the SKZCAM MRCC input strings for the adsorbate-slab complex, adsorbate, and slab in the same fashion as for ORCAInputGenerator.orcablocks
        self.skzcam_input_str: BlockInfo = {
            "adsorbate_slab": "",
            "adsorbate": "",
            "slab": "",
        }

        # Initialize the dictionary with keyword and values pairs for MRCC input
        self.skzcam_input_dict: MRCCInputDict = {
            "adsorbate_slab": {},
            "adsorbate": {},
            "slab": {},
        }

    def generate_input(self) -> MRCCInputDict:
        """
        Creates the mrccinput input for the MRCC ASE calculator.

        Returns
        -------
        MRCCInputDict
            A dictionary of key-value pairs (to be put in 'mrccinput' parameter) for the adsorbate-slab complex, the adsorbate, and the slab.
        """

        def _convert_input_str_to_dict(input_str: str) -> dict[str, str]:
            """
            Convert the SKZCAM input string to a dictionary.

            Parameters
            ----------
            input_str
                The SKZCAM input string containing all the input parameters for the SKZCAM protocol (i.e., basis, ecp, geometry, point charges)

            Returns
            -------
            dict[str,str]
                The SKZCAM input as a dictionary where each key is the input parameter and the value is the value of that parameter.
            """

            input_dict = {}

            key = None

            for line in input_str.split("\n"):
                if "=" in line:
                    key = line.split("=")[0]
                    input_dict[key] = line.split("=")[1]
                elif key is not None:
                    input_dict[key] += "\n" + line

            return input_dict

        # Create the blocks for the basis sets (basis, basis_sm, dfbasis_scf, dfbasis_cor, ecp)
        self._generate_basis_ecp_block()

        # Create the blocks for the coordinates
        self._generate_coords_block()

        # Create the point charge block and add it to the adsorbate-slab complex and slab blocks
        point_charge_block = self._generate_point_charge_block()
        self.skzcam_input_str["adsorbate_slab"] += point_charge_block
        self.skzcam_input_str["slab"] += point_charge_block

        # Convert the input string to a dictionary
        self.skzcam_input_dict["adsorbate_slab"] = _convert_input_str_to_dict(
            self.skzcam_input_str["adsorbate_slab"]
        )
        self.skzcam_input_dict["adsorbate"] = _convert_input_str_to_dict(
            self.skzcam_input_str["adsorbate"]
        )
        self.skzcam_input_dict["slab"] = _convert_input_str_to_dict(
            self.skzcam_input_str["slab"]
        )

        return self.skzcam_input_dict

    def _generate_basis_ecp_block(self) -> None:
        """
        Generates the basis and ECP block for the MRCC input file.

        Returns
        -------
        None
        """

        # Helper to generate basis strings for MRCC
        def _create_basis_block(quantum_region, ecp_region=None):

            atomtype_ecp = 'ecp=atomtype\n'
            for atom in quantum_region:
                atomtype_ecp += f"{self.element_info[atom.symbol]['ecp']}\n"
            if ecp_region is not None:
                atomtype_ecp += "capECP\n" * len(ecp_region)

            atomtype_ecp += '\n'

            return f"""
basis_sm=atomtype
{self._create_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: 'def2-SVP' for element in self.element_info})}

basis=atomtype
{self._create_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: self.element_info[element]['basis'] for element in self.element_info})}

dfbasis_scf=atomtype
{self._create_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: self.element_info[element]['ri_scf_basis'] for element in self.element_info})}

dfbasis_cor=atomtype
{self._create_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: self.element_info[element]['ri_cwft_basis'] for element in self.element_info})}

""" + atomtype_ecp

        if self.include_cp:
            self.skzcam_input_str["adsorbate_slab"] += _create_basis_block(
                quantum_region=self.adsorbate_slab_cluster, ecp_region=self.ecp_region
            )
            self.skzcam_input_str["slab"] += _create_basis_block(
                quantum_region=self.adsorbate_slab_cluster, ecp_region=self.ecp_region
            )
            self.skzcam_input_str["adsorbate"] += _create_basis_block(
                quantum_region=self.adsorbate_slab_cluster, ecp_region=None
            )
        else:
            self.skzcam_input_str["adsorbate_slab"] += _create_basis_block(
                quantum_region=self.adsorbate_slab_cluster, ecp_region=self.ecp_region
            )
            self.skzcam_input_str["slab"] += _create_basis_block(
                quantum_region=self.slab_cluster, ecp_region=self.ecp_region
            )
            self.skzcam_input_str["adsorbate"] += _create_basis_block(
                quantum_region=self.adsorbate_cluster, ecp_region=None
            )

    def _create_atomtype_basis(
        self,
        quantum_region: Atoms,
        element_basis_info: dict[ElementStr, str],
        ecp_region: Atoms | None = None,
    ) -> str:
        """
        Creates a column for the basis set for each atom in the Atoms object, given by element_info.

        Parameters
        ----------
        quantum_region
            The ASE Atoms object containing the atomic coordinates of the quantum cluster region (could be the adsorbate-slab complex, slab or adsorbate by itself).
        element_basis_info
            A dictionary with elements as keys which gives the basis set for each element.
        ecp_region
            The ASE atoms object containing the atomic coordinates of the capped ECP region.

        Returns
        -------
        str
            The basis set for each atom in the Atoms object given as a column (of size N, where N is the number of atoms).
        """

        basis_str = ""
        for atom in quantum_region:
            basis_str += f"{element_basis_info[atom.symbol]}\n"
        if ecp_region is not None:
            basis_str += "no-basis-set\n" * len(ecp_region)

        return basis_str

    def _generate_coords_block(self) -> None:
        """
        Generates the coordinates block for the MRCC input file. This includes the coordinates of the quantum cluster, the ECP region, and the point charges. It will return three strings for the adsorbate-slab complex, adsorbate and slab.

        Returns
        -------
        None
        """

        # Get the charge of the quantum cluster
        charge = int(sum(self.adsorbate_slab_cluster.get_array("oxi_states")))

        # Get the total number of core electrons for the quantum cluster
        core = {
            "adsorbate_slab": sum(
                [
                    self.element_info[atom.symbol]["core"]
                    for atom in self.adsorbate_slab_cluster
                ]
            ),
            "adsorbate": sum(
                [
                    self.element_info[atom.symbol]["core"]
                    for atom in self.adsorbate_cluster
                ]
            ),
            "slab": sum(
                [self.element_info[atom.symbol]["core"] for atom in self.slab_cluster]
            ),
        }

        # Add the charge and core electron information to skzcam_input_str
        self.skzcam_input_str["adsorbate_slab"] += f"""charge={charge}
mult={self.multiplicities['adsorbate_slab']}
core={int(core['adsorbate_slab']/2)}
unit=angs
geom=xyz
"""
        self.skzcam_input_str["adsorbate"] += f"""charge=0
mult={self.multiplicities['adsorbate']}
core={int(core['adsorbate']/2)}
unit=angs
geom=xyz
"""
        self.skzcam_input_str["slab"] += f"""charge={charge}
mult={self.multiplicities['slab']}
core={int(core['slab']/2)}
unit=angs
geom=xyz
"""
        # Create the atom coordinates block for the adsorbate-slab cluster, ECP region
        adsorbate_slab_coords_block = ""
        for atom in self.adsorbate_slab_cluster:
            adsorbate_slab_coords_block += create_atom_coord_string(atom=atom)

        ecp_region_block = ""
        for ecp_atom in self.ecp_region:
            ecp_region_block += create_atom_coord_string(atom=ecp_atom)

        # Set the number of atoms for each system. This would be the number of atoms in the quantum cluster plus the number of atoms in the ECP region. If include_cp is True, then the number of atoms in the quantum cluster is the number of atoms in the adsorbate-slab complex for both the adsorbate and slab.
        if self.include_cp:
            self.skzcam_input_str["adsorbate_slab"] += (
                f"{len(self.adsorbate_slab_cluster) + len(self.ecp_region)}\n\n"
            )
            self.skzcam_input_str["adsorbate_slab"] += (
                adsorbate_slab_coords_block + ecp_region_block
            )

            self.skzcam_input_str["adsorbate"] += (
                f"{len(self.adsorbate_slab_cluster)}\n\n"
            )
            self.skzcam_input_str["adsorbate"] += adsorbate_slab_coords_block

            self.skzcam_input_str["slab"] += (
                f"{len(self.adsorbate_slab_cluster) + len(self.ecp_region)}\n\n"
            )
            self.skzcam_input_str["slab"] += (
                adsorbate_slab_coords_block + ecp_region_block
            )

            for system in ["adsorbate_slab", "adsorbate", "slab"]:
                self.skzcam_input_str[system] += "\nghost=serialno\n"
                # Add the ghost atoms for the counterpoise correction in the adsorbate and slab
                if system == "adsorbate":
                    self.skzcam_input_str[system] += ",".join(
                        [str(atom_idx + 1) for atom_idx in self.slab_indices]
                    )
                elif system == "slab":
                    self.skzcam_input_str[system] += ",".join(
                        [str(atom_idx + 1) for atom_idx in self.adsorbate_indices]
                    )
                self.skzcam_input_str[system] += "\n\n"
        else:
            self.skzcam_input_str["adsorbate_slab"] += (
                f"{len(self.adsorbate_slab_cluster) + len(self.ecp_region)}\n\n"
            )
            self.skzcam_input_str["adsorbate_slab"] += (
                adsorbate_slab_coords_block + ecp_region_block
            )

            self.skzcam_input_str["adsorbate"] += f"{len(self.adsorbate_cluster)}\n\n"
            for atom in self.adsorbate_cluster:
                self.skzcam_input_str["adsorbate"] += create_atom_coord_string(
                    atom=atom
                )

            self.skzcam_input_str["slab"] += (
                f"{len(self.slab_cluster) + len(self.ecp_region)}\n\n"
            )
            for atom in self.slab_cluster:
                self.skzcam_input_str["slab"] += create_atom_coord_string(atom=atom)
            self.skzcam_input_str["slab"] += ecp_region_block

    def _generate_point_charge_block(self) -> str:
        """
        Create the point charge block for the MRCC input file. This requires the embedded_cluster Atoms object containing both atom_type and oxi_states arrays, as well as the indices of the quantum cluster and ECP region. Such arrays are created by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.

        Returns
        -------
        str
            The point charge block for the MRCC input file.
        """

        # Get the oxi_states arrays from the embedded_cluster
        oxi_states = self.adsorbate_slab_embedded_cluster.get_array("oxi_states")

        # Get the number of point charges for this system. There is a point charge associated with each capped ECP as well.
        pc_region_indices = [
            atom.index
            for atom in self.adsorbate_slab_embedded_cluster
            if atom.index not in self.quantum_cluster_indices
        ]

        num_pc = len(pc_region_indices)
        pc_block = f"qmmm=Amber\npointcharges\n{num_pc}\n"

        # Add the ecp_region indices
        for i in pc_region_indices:
            position = self.adsorbate_slab_embedded_cluster[i].position
            pc_block += f"  {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f} {oxi_states[i]:-16.11f}\n"

        return pc_block


class ORCAInputGenerator:
    """
    A class to generate the SKZCAM input for the ORCA ASE calculator.

    Attributes
    ----------
    adsorbate_slab_embedded_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
    quantum_cluster_indices
        A list containing the indices of the atoms in one quantum cluster. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
    ecp_region_indices
        A list containing the indices of the atoms in the ECP region of one quantum cluster. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    include_cp
        If True, the coords strings will include the counterpoise correction (i.e., ghost atoms) for the adsorbate and slab.
    multiplicities
        The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.
    pal_nprocs_block
        A dictionary with the number of processors for the PAL block as 'nprocs' and the maximum memory-per-core in megabytes blocks as 'maxcore'.
    method_block
        A dictionary that contains the method block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
    scf_block
        A dictionary that contains the SCF block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
    ecp_info
        A dictionary with the ECP data (in ORCA format) for the cations in the ECP region. The keys are the element symbols and the values are the ECP data.
    adsorbate_slab_cluster
        The ASE Atoms object for the quantum cluster of the adsorbate-slab complex.
    ecp_region
        The ASE Atoms object for the ECP region.
    adsorbate_indices
        The indices of the adsorbates from the adsorbate_slab_cluster quantum cluster.
    slab_indices
        The indices of the slab from the adsorbate_slab_cluster quantum cluster.
    adsorbate_cluster
        The ASE Atoms object for the quantum cluster of the adsorbate.
    slab_cluster
        The ASE Atoms object for the quantum cluster of the slab.
    orcablocks
        The ORCA input block (to be put in 'orcablocks' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.

    """

    def __init__(
        self,
        adsorbate_slab_embedded_cluster: Atoms,
        quantum_cluster_indices: list[int],
        ecp_region_indices: list[int],
        element_info: dict[ElementStr, ElementInfo] | None = None,
        include_cp: bool = True,
        multiplicities: MultiplicityDict | None = None,
        pal_nprocs_block: dict[str, int] | None = None,
        method_block: dict[str, str] | None = None,
        scf_block: dict[str, str] | None = None,
        ecp_info: dict[ElementStr, str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        adsorbate_slab_embedded_cluster
            The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        quantum_cluster_indices
            A list containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        ecp_region_indices
            A list containing the indices of the atoms in each ECP region. These indices are provided by the [quacc.atoms.skzcam.CreateSKZCAMClusters][] class.
        element_info
            A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
        include_cp
            If True, the coords strings will include the counterpoise correction (i.e., ghost atoms) for the adsorbate and slab.
        multiplicities
            The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.
        pal_nprocs_block
            A dictionary with the number of processors for the PAL block as 'nprocs' and the maximum memory-per-core in megabytes blocks as 'maxcore'.
        method_block
            A dictionary that contains the method block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
        scf_block
            A dictionary that contains the SCF block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
        ecp_info
            A dictionary with the ECP data (in ORCA format) for the cations in the ECP region.

        Returns
        -------
        None
        """

        self.adsorbate_slab_embedded_cluster = adsorbate_slab_embedded_cluster
        self.quantum_cluster_indices = quantum_cluster_indices
        self.ecp_region_indices = ecp_region_indices
        self.element_info = element_info
        self.include_cp = include_cp
        self.multiplicities = (
            {"adsorbate_slab": 1, "adsorbate": 1, "slab": 1}
            if multiplicities is None
            else multiplicities
        )
        self.pal_nprocs_block = pal_nprocs_block
        self.method_block = method_block
        self.scf_block = scf_block
        self.ecp_info = ecp_info

        # Check that none of the indices in quantum_cluster_indices are in ecp_region_indices
        if not np.all(
            [x not in self.ecp_region_indices for x in self.quantum_cluster_indices]
        ):
            raise ValueError(
                "An atom in the quantum cluster is also in the ECP region."
            )

        # Create the adsorbate-slab complex quantum cluster and ECP region cluster
        self.adsorbate_slab_cluster: Atoms = self.adsorbate_slab_embedded_cluster[
            self.quantum_cluster_indices
        ]
        self.ecp_region: Atoms = self.adsorbate_slab_embedded_cluster[
            self.ecp_region_indices
        ]

        # Get the indices of the adsorbates from the quantum cluster
        self.adsorbate_indices: list[int] = [
            i
            for i in range(len(self.adsorbate_slab_cluster))
            if self.adsorbate_slab_cluster.get_array("atom_type")[i] == "adsorbate"
        ]
        # Get the indices of the slab from the quantum cluster
        self.slab_indices: list[int] = [
            i
            for i in range(len(self.adsorbate_slab_cluster))
            if self.adsorbate_slab_cluster.get_array("atom_type")[i] != "adsorbate"
        ]

        # Create the adsorbate and slab quantum clusters
        self.adsorbate_cluster: Atoms = self.adsorbate_slab_cluster[
            self.adsorbate_indices
        ]
        self.slab_cluster: Atoms = self.adsorbate_slab_cluster[self.slab_indices]

        # Initialize the orcablocks input strings for the adsorbate-slab complex, adsorbate, and slab
        self.orcablocks: BlockInfo = {"adsorbate_slab": "", "adsorbate": "", "slab": ""}

    def generate_input(self) -> BlockInfo:
        """
        Creates the orcablocks input for the ORCA ASE calculator.

        Returns
        -------
        BlockInfo
            The ORCA input block (to be put in 'orcablocks' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
        """

        # First generate the preamble block
        self._generate_preamble_block()

        # Create the blocks for the coordinates
        self._generate_coords_block()

        # Combine the blocks
        return self.orcablocks

    def create_point_charge_file(self, pc_file: str | Path) -> None:
        """
        Create a point charge file that can be read by ORCA. This requires the embedded_cluster Atoms object containing both atom_type and oxi_states arrays, as well as the indices of the quantum cluster and ECP region.

        Parameters
        ----------
        pc_file
            A file containing the point charges to be written by ORCA.

        Returns
        -------
        None
        """

        # Get the oxi_states arrays from the embedded_cluster
        oxi_states = self.adsorbate_slab_embedded_cluster.get_array("oxi_states")

        # Get the number of point charges for this system
        total_indices = self.quantum_cluster_indices + self.ecp_region_indices
        num_pc = len(self.adsorbate_slab_embedded_cluster) - len(total_indices)
        counter = 0
        with Path.open(pc_file, "w") as f:
            # Write the number of point charges first
            f.write(f"{num_pc}\n")
            for i in range(len(self.adsorbate_slab_embedded_cluster)):
                if i not in total_indices:
                    counter += 1
                    position = self.adsorbate_slab_embedded_cluster[i].position
                    if counter != num_pc:
                        f.write(
                            f"{oxi_states[i]:-16.11f} {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f}\n"
                        )
                    else:
                        f.write(
                            f"{oxi_states[i]:-16.11f} {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f}"
                        )

    def _generate_coords_block(self) -> None:
        """
        Generates the coordinates block for the ORCA input file. This includes the coordinates of the quantum cluster, the ECP region, and the point charges. It will return three strings for the adsorbate-slab complex, adsorbate and slab.

        Returns
        -------
        None
        """

        # Get the charge of the adsorbate_slab cluster
        charge = int(sum(self.adsorbate_slab_cluster.get_array("oxi_states")))

        # Add the coords strings for the adsorbate-slab complex, adsorbate, and slab
        self.orcablocks["adsorbate_slab"] += f"""%coords
CTyp xyz
Mult {self.multiplicities['adsorbate_slab']}
Units angs
Charge {charge}
coords
"""
        self.orcablocks["adsorbate"] += f"""%coords
CTyp xyz
Mult {self.multiplicities['adsorbate']}
Units angs
Charge 0
coords
"""
        self.orcablocks["slab"] += f"""%coords
CTyp xyz
Mult {self.multiplicities['slab']}
Units angs
Charge {charge}
coords
"""

        for i, atom in enumerate(self.adsorbate_slab_cluster):
            # Create the coords section for the adsorbate-slab complex
            self.orcablocks["adsorbate_slab"] += create_atom_coord_string(atom=atom)

            # Create the coords section for the adsorbate and slab
            if i in self.adsorbate_indices:
                self.orcablocks["adsorbate"] += create_atom_coord_string(atom=atom)
                if self.include_cp:
                    self.orcablocks["slab"] += create_atom_coord_string(
                        atom=atom, is_ghost_atom=True
                    )
            elif i in self.slab_indices:
                self.orcablocks["slab"] += create_atom_coord_string(atom=atom)
                if self.include_cp:
                    self.orcablocks["adsorbate"] += create_atom_coord_string(
                        atom=atom, is_ghost_atom=True
                    )

        # Create the coords section for the ECP region
        ecp_region_coords_section = ""
        for i, atom in enumerate(self.ecp_region):
            atom_ecp_info = self._format_ecp_info(
                atom_ecp_info=self.ecp_info[atom.symbol]
            )
            ecp_region_coords_section += create_atom_coord_string(
                atom=atom,
                atom_ecp_info=atom_ecp_info,
                pc_charge=self.ecp_region.get_array("oxi_states")[i],
            )

        # Add the ECP region coords section to the ads_slab_coords string
        self.orcablocks["adsorbate_slab"] += f"{ecp_region_coords_section}end\nend\n"
        self.orcablocks["slab"] += f"{ecp_region_coords_section}end\nend\n"
        self.orcablocks["adsorbate"] += "end\nend\n"

    def _format_ecp_info(self, atom_ecp_info: str) -> str:
        """
        Formats the ECP info so that it can be inputted to ORCA without problems.

        Parameters
        ----------
        atom_ecp_info
            The ECP info for a single atom.

        Returns
        -------
        str
            The formatted ECP info.
        """
        # Find the starting position of "NewECP" and "end"
        start_pos = atom_ecp_info.lower().find("newecp")
        end_pos = atom_ecp_info.lower().find("end", start_pos)

        start_pos += len("NewECP")

        # If "NewECP" or "end" is not found, then we assume that ecp_info has been given without these lines but in the correct format
        if start_pos == -1 or end_pos == -1:
            raise ValueError("ECP info does not contain 'NewECP' or 'end' keyword.")

        # Extract content between "NewECP" and "end", exclusive of "end", then add correctly formatted "NewECP" and "end"
        return f"NewECP\n{atom_ecp_info[start_pos:end_pos].strip()}\nend\n"

    def _generate_preamble_block(self) -> str:
        """
        From the quantum cluster Atoms object, generate the ORCA input preamble for the basis, method, pal, and scf blocks.

        Returns
        -------
        None
        """

        # Get the set of element symbols from the quantum cluster
        element_symbols = list(set(self.adsorbate_slab_cluster.get_chemical_symbols()))
        element_symbols.sort()

        # Check all element symbols are provided in element_info keys
        if self.element_info is not None and not all(
            element in self.element_info for element in element_symbols
        ):
            raise ValueError(
                "Not all element symbols are provided in the element_info dictionary."
            )

        # Initialize preamble_input
        preamble_input = ""

        # Add the pal_nprocs_block
        if self.pal_nprocs_block is not None:
            preamble_input += f"%pal nprocs {self.pal_nprocs_block['nprocs']} end\n"
            preamble_input += f"%maxcore {self.pal_nprocs_block['maxcore']} end\n"

        # Add pointcharge file to read. It will be assumed that it is in the same folder as the input file
        preamble_input += '%pointcharges "orca.pc"\n'

        # Make the method block
        if self.method_block is not None and self.element_info is not None:
            preamble_input += "%method\n"
        # Iterate through the keys of method_block and add key value
        if self.method_block is not None:
            for key in self.method_block:
                preamble_input += f"{key} {self.method_block[key]}\n"
        # Iterate over the core value for each element (if it has been given)
        if self.element_info is not None:
            for element in element_symbols:
                if "core" in self.element_info[element]:
                    preamble_input += (
                        f"NewNCore {element} {self.element_info[element]['core']} end\n"
                    )
        if self.method_block is not None and self.element_info is not None:
            preamble_input += "end\n"

        # Make the basis block

        # First check if the basis key is the same for all elements. We use """ here because an option for these keys is "AutoAux"
        if self.element_info is not None:
            preamble_input += "%basis\n"
            if (
                len(
                    {self.element_info[element]["basis"] for element in element_symbols}
                )
                == 1
            ):
                preamble_input += (
                    f"""Basis {self.element_info[element_symbols[0]]['basis']}\n"""
                )
            else:
                for element in element_symbols:
                    element_basis = self.element_info[element]["basis"]
                    preamble_input += f"""NewGTO {element} "{element_basis}" end\n"""

            # Do the same for ri_scf_basis and ri_cwft_basis.
            if (
                len(
                    {
                        self.element_info[element]["ri_scf_basis"]
                        for element in element_symbols
                    }
                )
                == 1
            ):
                preamble_input += (
                    f"""Aux {self.element_info[element_symbols[0]]['ri_scf_basis']}\n"""
                )
            else:
                for element in element_symbols:
                    element_basis = self.element_info[element]["ri_scf_basis"]
                    preamble_input += f'NewAuxJGTO {element} "{element_basis}" end\n'

            if (
                len(
                    list(
                        {
                            self.element_info[element]["ri_cwft_basis"]
                            for element in element_symbols
                        }
                    )
                )
                == 1
            ):
                preamble_input += f"""AuxC {self.element_info[element_symbols[0]]['ri_cwft_basis']}\n"""
            else:
                for element in element_symbols:
                    element_basis = self.element_info[element]["ri_cwft_basis"]
                    preamble_input += (
                        f"""NewAuxCGTO {element} "{element_basis}" end\n"""
                    )

            preamble_input += "end\n"

        # Write the scf block
        if self.scf_block is not None:
            preamble_input += "%scf\n"
            for key in self.scf_block:
                preamble_input += f"""{key} {self.scf_block[key]}\n"""
            preamble_input += "end\n"

        # Add preamble_input to the orcablocks for the adsorbate-slab complex, adsorbate, and slab
        self.orcablocks["adsorbate_slab"] += preamble_input
        self.orcablocks["adsorbate"] += "\n".join([line for line in preamble_input.splitlines() if "orca.pc" not in line])
        self.orcablocks["slab"] += preamble_input


class CreateSKZCAMClusters:
    """
    A class to create the quantum clusters and ECP regions for the SKZCAM protocol.

    Attributes
    ----------
    adsorbate_indices
        The indices of the atoms that make up the adsorbate molecule.
    slab_center_indices
        The indices of the atoms that make up the 'center' of the slab right beneath the adsorbate.
    slab_indices
        The indices of the atoms that make up the slab.
    atom_oxi_states
        A dictionary with the element symbol as the key and its oxidation state as the value.
    adsorbate_slab_file
        The path to the file containing the adsorbate molecule on the surface slab. It can be in any format that ASE can read.
    pun_file
        The path to the .pun file containing the atomic coordinates and charges of the adsorbate-slab complex. This file should be generated by ChemShell. If it is None, then ChemShell wil be used to create this file.
    adsorbate
        The ASE Atoms object containing the atomic coordinates of the adsorbate.
    slab
        The ASE Atoms object containing the atomic coordinates of the slab.
    adsorbate_slab
        The ASE Atoms object containing the atomic coordinates of the adsorbate-slab complex.
    adsorbate_slab_embedded_cluster
        The ASE Atoms object containing the atomic coordinates, atomic charges and atom type (i.e., point charge or cation/anion) from the .pun file for the embedded cluster of the adsorbate-slab complex.
    slab_embedded_cluster
        The ASE Atoms object containing the atomic coordinates, atomic charges and atom type (i.e., point charge or cation/anion) from the .pun file for the embedded cluster of the slab.
    quantum_cluster_indices_set
        A list of lists of indices of the atoms in the set of quantum clusters created by the SKZCAM protocol
    ecp_region_indices_set
        A list of lists of indices of the atoms in the ECP region for the set of quantum clusters created by the SKZCAM protocol
    """

    def __init__(
        self,
        adsorbate_indices: list[int],
        slab_center_indices: list[int],
        atom_oxi_states: dict[str, int],
        adsorbate_slab_file: str | Path | None = None,
        pun_file: str | Path | None = None,
    ) -> None:
        """
        Parameters
        ----------
        adsorbate_indices
            The indices of the atoms that make up the adsorbate molecule.
        slab_center_indices
            The indices of the atoms that make up the 'center' of the slab right beneath the adsorbate.
        atom_oxi_states
            A dictionary with the element symbol as the key and its oxidation state as the value.
        adsorbate_slab_file
            The path to the file containing the adsorbate molecule on the surface slab. It can be in any format that ASE can read.
        pun_file
            The path to the .pun file containing the atomic coordinates and charges of the adsorbate-slab complex. This file should be generated by ChemShell. If it is None, then ChemShell wil be used to create this file.

        Returns
        -------
        None
        """

        self.adsorbate_indices = adsorbate_indices
        self.slab_center_indices = slab_center_indices
        self.slab_indices = None  # This will be set later
        self.atom_oxi_states = atom_oxi_states
        self.adsorbate_slab_file = adsorbate_slab_file
        self.pun_file = pun_file

        # Check that the adsorbate_indices and slab_center_indices are not the same
        if any(x in self.adsorbate_indices for x in self.slab_center_indices):
            raise ValueError(
                "The adsorbate and slab center indices cannot be the same."
            )

        # Check that the adsorbate_slab_file and pun_file are not both None
        if self.adsorbate_slab_file is None and self.pun_file is None:
            raise ValueError(
                "Either the adsorbate_slab_file or pun_file must be provided."
            )

        # Initialize the adsorbate, slab and adsorbate_slab Atoms object which contains the adsorbate, slab and adsorbate-slab complex respectively
        self.adsorbate: Atoms | None
        self.slab: Atoms | None
        self.adsorbate_slab: Atoms | None

        # Initialize the embedded_adsorbate_slab_cluster, and embedded_slab_cluster Atoms object which are the embedded cluster for the adsorbate-slab complex and slab respectively
        self.adsorbate_slab_embedded_cluster: Atoms | None = None
        self.slab_embedded_cluster: Atoms | None = None

        # Initialize the quantum cluster indices and ECP region indices
        self.quantum_cluster_indices_set: list[list[int]] | None = None
        self.ecp_region_indices_set: list[list[int]] | None = None

    def convert_slab_to_atoms(self) -> None:
        """
        Read the file containing the periodic slab and adsorbate (geometry optimized) and format the resulting Atoms object to be used to create a .pun file in ChemShell.

        Returns
        -------
        None
        """

        # Get the necessary information for the cluster from a provided slab file (in any format that ASE can read)
        adsorbate_slab = read(self.adsorbate_slab_file)

        # Find indices (within adsorbate_slab) of the slab
        slab_indices = self.slab_center_indices + [
            i
            for i, _ in enumerate(adsorbate_slab)
            if i not in (self.adsorbate_indices + self.slab_center_indices)
        ]

        # Create adsorbate and slab from adsorbate_slab
        slab = adsorbate_slab[slab_indices]
        adsorbate = adsorbate_slab[self.adsorbate_indices]

        adsorbate.translate(-slab[0].position)
        slab.translate(-slab[0].position)

        # Get the relative distance of the adsorbate from the first center atom of the slab as defined in the slab_center_indices
        adsorbate_vector_from_slab = adsorbate[0].position - slab[0].position

        # Get the center of the cluster from the slab_center_indices
        slab_center_position = slab[
            : len(self.slab_center_indices)
        ].get_positions().sum(axis=0) / len(self.slab_center_indices)

        # Add the height of the adsorbate from the slab along the z-direction relative to the slab_center
        adsorbate_com_z_disp = (
            adsorbate.get_center_of_mass()[2] - slab_center_position[2]
        )

        center_position = (
            np.array([0.0, 0.0, adsorbate_com_z_disp]) + slab_center_position
        )

        self.adsorbate = adsorbate
        self.slab = slab
        self.adsorbate_slab = adsorbate_slab
        self.adsorbate_vector_from_slab = adsorbate_vector_from_slab
        self.center_position = center_position

    @requires(has_chemshell, "ChemShell is not installed")
    def run_chemshell(
        self,
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
        filepath
            The location where the ChemShell output files will be written.
        chemsh_radius_active
            The radius of the active region in Angstroms. This 'active' region is simply region where the charge fitting is performed to ensure correct Madelung potential; it can be a relatively large value.
        chemsh_radius_cluster
            The radius of the total embedded cluster in Angstroms.
        chemsh_bq_layer
            The height above the surface to place some additional fitting point charges in Angstroms; simply for better reproduction of the electrostatic potential close to the adsorbate.
        write_xyz_file
            Whether to write an XYZ file of the cluster for visualisation.

        Returns
        -------
        None
        """
        from chemsh.io.tools import convert_atoms_to_frag

        # Convert ASE Atoms to ChemShell Fragment object
        slab_frag = convert_atoms_to_frag(self.slab, connect_mode="ionic", dim="2D")

        # Add the atomic charges to the fragment
        slab_frag.addCharges(self.atom_oxi_states)

        # Create the chemshell cluster (i.e., add electrostatic fitting charges) from the fragment
        chemsh_slab_embedded_cluster = slab_frag.construct_cluster(
            origin=0,
            radius_cluster=chemsh_radius_cluster / Bohr,
            radius_active=chemsh_radius_active / Bohr,
            bq_layer=chemsh_bq_layer / Bohr,
            adjust_charge="coordination_scaled",
        )

        # Save the final cluster to a .pun file
        chemsh_slab_embedded_cluster.save(
            filename=Path(filepath).with_suffix(".pun"), fmt="pun"
        )
        self.pun_file = Path(filepath).with_suffix(".pun")

        if write_xyz_file:
            # XYZ for visualisation
            chemsh_slab_embedded_cluster.save(
                filename=Path(filepath).with_suffix(".xyz"), fmt="xyz"
            )

    def run_skzcam(
        self,
        shell_max: int = 10,
        shell_width: float = 0.1,
        bond_dist: float = 2.5,
        ecp_dist: float = 6.0,
        write_clusters: bool = False,
        write_clusters_path: str | Path = ".",
        write_include_ecp: bool = False,
    ) -> SKZCAMOutput:
        """
        From a provided .pun file (generated by ChemShell), this function creates quantum clusters using the SKZCAM protocol. It will return the embedded cluster Atoms object and the indices of the atoms in the quantum clusters and the ECP region. The number of clusters created is controlled by the rdf_max parameter.

        Parameters
        ----------
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
        dict[str, Atoms]
            A dictionary containing the embedded cluster Atoms object of the adsorbate-slab complex, accessed by the key 'adsorbate_slab_embedded_cluster', and the indices of the atoms in the quantum clusters and the ECP region, accessed by the keys 'quantum_cluster_indices_set' and 'ecp_region_indices_set' respectively.
        """

        # Read the .pun file and create the embedded_cluster Atoms object
        self.slab_embedded_cluster = self._convert_pun_to_atoms(pun_file=self.pun_file)

        # Get distances of all atoms from the cluster center
        atom_center_distances = _get_atom_distances(
            atoms=self.slab_embedded_cluster, center_position=self.center_position
        )

        # Determine the cation shells from the center of the embedded cluster
        _, cation_shells_idx = self._find_cation_shells(
            slab_embedded_cluster=self.slab_embedded_cluster,
            distances=atom_center_distances,
            shell_width=shell_width,
        )

        # Create the distance matrix for the embedded cluster
        slab_embedded_cluster_all_dist = self.slab_embedded_cluster.get_all_distances()

        # Create the anion coordination list for each cation shell
        anion_coord_idx = []
        for shell_idx in range(shell_max):
            shell_indices = cation_shells_idx[shell_idx]
            anion_coord_idx += [
                self._get_anion_coordination(
                    slab_embedded_cluster=self.slab_embedded_cluster,
                    cation_shell_indices=shell_indices,
                    dist_matrix=slab_embedded_cluster_all_dist,
                    bond_dist=bond_dist,
                )
            ]

        # Create the quantum clusters by summing up the indices of the cations and their coordinating anions
        slab_quantum_cluster_indices_set = []
        dummy_cation_indices = []
        dummy_anion_indices = []
        for shell_idx in range(shell_max):
            dummy_cation_indices += cation_shells_idx[shell_idx]
            dummy_anion_indices += anion_coord_idx[shell_idx]
            slab_quantum_cluster_indices_set += [
                list(set(dummy_cation_indices + dummy_anion_indices))
            ]

        # Get the ECP region for each quantum cluster
        slab_ecp_region_indices_set = self._get_ecp_region(
            slab_embedded_cluster=self.slab_embedded_cluster,
            quantum_cluster_indices_set=slab_quantum_cluster_indices_set,
            dist_matrix=slab_embedded_cluster_all_dist,
            ecp_dist=ecp_dist,
        )

        # Create the adsorbate_slab_embedded_cluster from slab_embedded_cluster and adsorbate atoms objects. This also sets the final quantum_cluster_indices_set and ecp_region_indices_set for the adsorbate_slab_embedded_cluster
        self._create_adsorbate_slab_embedded_cluster(
            quantum_cluster_indices_set=slab_quantum_cluster_indices_set,
            ecp_region_indices_set=slab_ecp_region_indices_set,
        )

        # Write the quantum clusters to files
        if write_clusters:
            for idx in range(len(self.quantum_cluster_indices_set)):
                quantum_atoms = self.adsorbate_slab_embedded_cluster[
                    self.quantum_cluster_indices_set[idx]
                ]
                if write_include_ecp:
                    ecp_atoms = self.adsorbate_slab_embedded_cluster[
                        self.ecp_region_indices_set[idx]
                    ]
                    ecp_atoms.set_chemical_symbols(np.array(["U"] * len(ecp_atoms)))
                    cluster_atoms = quantum_atoms + ecp_atoms
                else:
                    cluster_atoms = quantum_atoms
                write(
                    Path(write_clusters_path, f"SKZCAM_cluster_{idx}.xyz"),
                    cluster_atoms,
                )

        return {
            "adsorbate_slab_embedded_cluster": self.adsorbate_slab_embedded_cluster,
            "quantum_cluster_indices_set": self.quantum_cluster_indices_set,
            "ecp_region_indices_set": self.ecp_region_indices_set,
        }

    def _convert_pun_to_atoms(self, pun_file: str | Path) -> Atoms:
        """
        Reads a .pun file and returns an ASE Atoms object containing the atomic coordinates,
        point charges/oxidation states, and atom types.

        Parameters
        ----------
        pun_file
            The path to the .pun file created by ChemShell to be read.

        Returns
        -------
        Atoms
            The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file.
            The `oxi_states` array contains the atomic charges, and the `atom_type` array contains the
            atom types (cation, anion, neutral).
        """

        # Create a dictionary containing the atom types and whether they are cations or anions
        atom_type_dict = {
            atom: "cation" if oxi_state > 0 else "anion" if oxi_state < 0 else "neutral"
            for atom, oxi_state in self.atom_oxi_states.items()
        }

        # Load the pun file as a list of strings
        with zopen(zpath(str(Path(pun_file)))) as f:
            raw_pun_file = [
                line.rstrip().decode("utf-8")
                if isinstance(line, bytes)
                else line.rstrip()
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
        for _, line in enumerate(raw_atom_positions):
            line_info = line.split()

            # Add the atom type to the atom_type_list
            if line_info[0] in atom_type_dict:
                atom_types.append(atom_type_dict[line_info[0]])
            elif line_info[0] == "F":
                atom_types.append("pc")
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

        slab_embedded_cluster = Atoms(numbers=atom_numbers, positions=atom_positions)

        # Center the embedded cluster so that atom index 0 is at the [0, 0, 0] position
        slab_embedded_cluster.translate(-slab_embedded_cluster[0].position)

        # Add the `oxi_states` and `atom_type` arrays to the Atoms object
        slab_embedded_cluster.set_array("oxi_states", np.array(charges))
        slab_embedded_cluster.set_array("atom_type", np.array(atom_types))

        return slab_embedded_cluster

    def _create_adsorbate_slab_embedded_cluster(
        self,
        quantum_cluster_indices_set: list[list[int]] | None = None,
        ecp_region_indices_set: list[list[int]] | None = None,
    ) -> None:
        """
        Insert the adsorbate into the embedded cluster and update the quantum cluster and ECP region indices.

        Parameters
        ----------
        quantum_cluster_indices_set
            A list of lists containing the indices of the atoms in each quantum cluster.
        ecp_region_indices_set
            A list of lists containing the indices of the atoms in the ECP region for each quantum cluster.

        Returns
        -------
        None
        """

        # Remove PBC from the adsorbate
        self.adsorbate.set_pbc(False)

        # Translate the adsorbate to the correct position relative to the slab
        self.adsorbate.translate(
            self.slab_embedded_cluster[0].position
            - self.adsorbate[0].position
            + self.adsorbate_vector_from_slab
        )

        # Set oxi_state and atom_type arrays for the adsorbate
        self.adsorbate.set_array("oxi_states", np.array([0.0] * len(self.adsorbate)))
        self.adsorbate.set_array(
            "atom_type", np.array(["adsorbate"] * len(self.adsorbate))
        )

        # Add the adsorbate to the embedded cluster
        self.adsorbate_slab_embedded_cluster = (
            self.adsorbate + self.slab_embedded_cluster
        )

        # Update the quantum cluster and ECP region indices
        if quantum_cluster_indices_set is not None:
            quantum_cluster_indices_set = [
                list(range(len(self.adsorbate)))
                + [idx + len(self.adsorbate) for idx in cluster]
                for cluster in quantum_cluster_indices_set
            ]
        if ecp_region_indices_set is not None:
            ecp_region_indices_set = [
                [idx + len(self.adsorbate) for idx in cluster]
                for cluster in ecp_region_indices_set
            ]

        self.quantum_cluster_indices_set = quantum_cluster_indices_set
        self.ecp_region_indices_set = ecp_region_indices_set

    def _find_cation_shells(
        self, slab_embedded_cluster: Atoms, distances: NDArray, shell_width: float = 0.1
    ) -> tuple[list[list[int]], list[list[int]]]:
        """
        Returns a list of lists containing the indices of the cations in each shell, based on distance from the embedded cluster center.
        This is achieved by clustering the data based on the DBSCAN clustering algorithm.

        Parameters
        ----------
        slab_embedded_cluster
            The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
        distances
            The distance of atoms from the cluster center.
        shell_width
            Defines the distance between atoms within shells; this is the maximum distance between any two atoms within the shell

        Returns
        -------
        list[list[int]]
            A list of lists containing the distance of the cation in each shell from the adsorbate.
        list[list[int]]
            A list of lists containing the indices of the cations in each shell.
        """

        # Define the empty list to store the cation shells
        shells_distances = []
        shells_indices = []

        # Sort the points by distance from the cluster center for the cations only
        distances_sorted = []
        distances_sorted_indices = []
        for i in np.argsort(distances):
            if slab_embedded_cluster.get_array("atom_type")[i] == "cation":
                distances_sorted.append(distances[i])
                distances_sorted_indices.append(i)

        current_point = distances_sorted[0]
        current_shell = [current_point]
        current_shell_idx = [distances_sorted_indices[0]]

        for idx, point in enumerate(distances_sorted[1:]):
            if point <= current_point + shell_width:
                current_shell.append(point)
                current_shell_idx.append(distances_sorted_indices[idx + 1])
            else:
                shells_distances.append(current_shell)
                shells_indices.append(current_shell_idx)
                current_shell = [point]
                current_shell_idx = [distances_sorted_indices[idx + 1]]
            current_point = point
        shells_distances.append(current_shell)
        shells_indices.append(current_shell_idx)

        return shells_distances, shells_indices

    def _get_anion_coordination(
        self,
        slab_embedded_cluster: Atoms,
        cation_shell_indices: list[int],
        dist_matrix: NDArray,
        bond_dist: float = 2.5,
    ) -> list[int]:
        """
        Returns a list of lists containing the indices of the anions coordinating the cation indices provided.

        Parameters
        ----------
        slab_embedded_cluster
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
                    and slab_embedded_cluster.get_array("atom_type")[idx] == "anion"
                )
            ]

        return list(set(anion_coord_indices))

    def _get_ecp_region(
        self,
        slab_embedded_cluster: Atoms,
        quantum_cluster_indices_set: list[int],
        dist_matrix: NDArray,
        ecp_dist: float = 6.0,
    ) -> list[list[int]]:
        """
        Returns a list of lists containing the indices of the atoms in the ECP region of the embedded cluster for each quantum cluster

        Parameters
        ----------
        slab_embedded_cluster
            The ASE Atoms object containing the atomic coordinates AND the atom types (i.e. cation or anion).
        quantum_cluster_indices_set
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

        ecp_region_indices_set = []
        dummy_cation_indices = []

        # Iterate over the quantum clusters and find the atoms within the ECP distance of each quantum cluster
        for cluster in quantum_cluster_indices_set:
            dummy_cation_indices += cluster
            cluster_ecp_region_idx = []
            for atom_idx in dummy_cation_indices:
                for idx, dist in enumerate(dist_matrix[atom_idx]):
                    # Check if the atom is within the ecp_dist region and is not in the quantum cluster and is a cation
                    if (
                        dist < ecp_dist
                        and idx not in dummy_cation_indices
                        and slab_embedded_cluster.get_array("atom_type")[idx]
                        == "cation"
                    ):
                        cluster_ecp_region_idx += [idx]

            ecp_region_indices_set += [list(set(cluster_ecp_region_idx))]

        return ecp_region_indices_set


def _get_atom_distances(atoms: Atoms, center_position: NDArray) -> NDArray:
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

    return np.array([np.linalg.norm(atom.position - center_position) for atom in atoms])


def create_atom_coord_string(
    atom: Atom,
    is_ghost_atom: bool = False,
    atom_ecp_info: str | None = None,
    pc_charge: float | None = None,
) -> str:
    """
    Creates a string containing the Atom symbol and coordinates for both MRCC and ORCA, with additional information for atoms in the ECP region as well as ghost atoms.

    Parameters
    ----------
    atom
        The ASE Atom (not Atoms) object containing the atomic coordinates.
    is_ghost_atom
        If True, then the atom is a ghost atom.
    atom_ecp_info
        If not None, then assume this is an atom in the ECP region and adds the ECP info.
    pc_charge
        The point charge value for the ECP region atom.

    Returns
    -------
    str
        The atom symbol and coordinates in the ORCA input file format.
    """

    # If ecp_info is not None and ghost_atom is True, raise an error
    if atom_ecp_info and is_ghost_atom:
        raise ValueError("ECP info cannot be provided for ghost atoms.")

    # Check that pc_charge is a float if atom_ecp_info is not None
    if atom_ecp_info and pc_charge is None:
        raise ValueError("Point charge value must be given for atoms with ECP info.")

    if is_ghost_atom:
        atom_coord_str = f"{(atom.symbol + ':').ljust(3)} {' '*16} {atom.position[0]:-16.11f} {atom.position[1]:-16.11f} {atom.position[2]:-16.11f}\n"
    elif atom_ecp_info is not None:
        atom_coord_str = f"{(atom.symbol + '>').ljust(3)} {pc_charge:-16.11f} {atom.position[0]:-16.11f} {atom.position[1]:-16.11f} {atom.position[2]:-16.11f}\n{atom_ecp_info}"
    else:
        atom_coord_str = f"{atom.symbol.ljust(3)} {' '*16} {atom.position[0]:-16.11f} {atom.position[1]:-16.11f} {atom.position[2]:-16.11f}\n"

    return atom_coord_str
