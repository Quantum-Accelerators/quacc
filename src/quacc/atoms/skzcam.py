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
    from ase.atom import Atom
    from numpy.typing import NDArray

    from quacc.types import BlockInfo, ElementInfo, ElementStr, MultiplicityDict


has_chemshell = find_spec("chemsh") is not None


def create_mrcc_eint_blocks(
    embedded_adsorbed_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
    element_info: dict[ElementStr, ElementInfo] | None = None,
    include_cp: bool = True,
    multiplicities: MultiplicityDict | None = None,
) -> BlockInfo:
    """
    Creates the orcablocks input for the MRCC ASE calculator.

    Parameters
    ----------
    embedded_adsorbed_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    quantum_cluster_indices
        A list containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list containing the indices of the atoms in each ECP region. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    include_cp
        If True, the coords strings will include the counterpoise correction for the adsorbate and slab.
    multiplicities
        The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.

    Returns
    -------
    BlockInfo
        The ORCA input block (to be put in 'orcablocks' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
    """

    # Create the blocks for the basis sets (basis, basis_sm, dfbasis_scf, dfbasis_cor, ecp)
    basis_ecp_block = generate_mrcc_basis_ecp_block(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=quantum_cluster_indices,
        ecp_region_indices=ecp_region_indices,
        element_info=element_info,
        include_cp=include_cp,
    )

    # Create the blocks for the coordinates
    coords_block = generate_mrcc_coords_block(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=quantum_cluster_indices,
        ecp_region_indices=ecp_region_indices,
        element_info=element_info,
        include_cp=include_cp,
        multiplicities=multiplicities,
    )

    # Create the point charge block
    point_charge_block = generate_mrcc_point_charge_block(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=quantum_cluster_indices,
        ecp_region_indices=ecp_region_indices,
    )

    # Combine the blocks
    return {
        "adsorbate_slab": basis_ecp_block["adsorbate_slab"]
        + coords_block["adsorbate_slab"]
        + point_charge_block,
        "adsorbate": basis_ecp_block["adsorbate"] + coords_block["adsorbate"],
        "slab": basis_ecp_block["slab"] + coords_block["slab"] + point_charge_block,
    }


def generate_mrcc_basis_ecp_block(
    embedded_adsorbed_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
    element_info: dict[ElementStr, ElementInfo] | None = None,
    include_cp: bool = True,
) -> BlockInfo:
    """
    Generates the basis and ECP block for the MRCC input file.

    Parameters
    ----------
    embedded_adsorbed_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file. This object is created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function and should contain the atom_types for the adsorbate and slab.
    quantum_cluster_indices
        A list of lists containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    include_cp
        If True, the coords strings will include the counterpoise correction (i.e., ghost atoms) for the adsorbate and slab.

    Returns
    -------
    BlockInfo
        The basis and ECP block for the MRCC input file for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
    """

    # Create the quantum cluster and ECP region cluster
    adsorbate_slab_cluster = embedded_adsorbed_cluster[quantum_cluster_indices]
    ecp_region = embedded_adsorbed_cluster[ecp_region_indices]

    # Get the indices of the adsorbates from the quantum cluster
    adsorbate_indices = [
        i
        for i in range(len(adsorbate_slab_cluster))
        if adsorbate_slab_cluster.get_array("atom_type")[i] == "adsorbate"
    ]

    # Get the indices of the slab from the quantum cluster
    slab_indices = [
        i
        for i in range(len(adsorbate_slab_cluster))
        if adsorbate_slab_cluster.get_array("atom_type")[i] != "adsorbate"
    ]

    adsorbate_cluster = adsorbate_slab_cluster[adsorbate_indices]
    slab_cluster = adsorbate_slab_cluster[slab_indices]

    # Helper to generate basis strings for MRCC
    def _create_basis_block(quantum_region, ecp_region=None):
        return f"""
basis_sm=atomtype
{create_mrcc_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: 'def2-SVP' for element in element_info})}

basis=atomtype
{create_mrcc_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: element_info[element]['basis'] for element in element_info})}

dfbasis_scf=atomtype
{create_mrcc_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: element_info[element]['ri_scf_basis'] for element in element_info})}

dfbasis_cor=atomtype
{create_mrcc_atomtype_basis(quantum_region=quantum_region, ecp_region=ecp_region, element_basis_info={element: element_info[element]['ri_cwft_basis'] for element in element_info})}
"""

    if include_cp:
        return {
            "adsorbate_slab": _create_basis_block(
                quantum_region=adsorbate_slab_cluster, ecp_region=ecp_region
            ),
            "slab": _create_basis_block(
                quantum_region=adsorbate_slab_cluster, ecp_region=ecp_region
            ),
            "adsorbate": _create_basis_block(
                quantum_region=adsorbate_slab_cluster, ecp_region=None
            ),
        }
    else:
        return {
            "adsorbate_slab": _create_basis_block(
                quantum_region=adsorbate_slab_cluster, ecp_region=ecp_region
            ),
            "slab": _create_basis_block(
                quantum_region=slab_cluster, ecp_region=ecp_region
            ),
            "adsorbate": _create_basis_block(
                quantum_region=adsorbate_cluster, ecp_region=None
            ),
        }


def create_mrcc_atomtype_basis(
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


def generate_mrcc_coords_block(
    embedded_adsorbed_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
    element_info: ElementInfo,
    include_cp: bool = True,
    multiplicities: MultiplicityDict | None = None,
) -> BlockInfo:
    """
    Generates the coordinates block for the MRCC input file. This includes the coordinates of the quantum cluster, the ECP region, and the point charges. It will return three strings for the adsorbate-slab complex, adsorbate and slab.

    Parameters
    ----------
    embedded_adsorbed_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This Atoms object is typically produced from [quacc.atoms.skzcam.create_skzcam_clusters][].
    quantum_cluster_indices
        A list containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list containing the indices of the atoms in each ECP region. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    include_cp
        If True, the coords strings will include the counterpoise correction for the adsorbate and slab.
    multiplicities
        The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.

    Returns
    -------
    BlockInfo
        The coordinates block for the ORCA input block (to be put in 'orcablocks' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
    """

    # Create the quantum cluster and ECP region cluster
    if multiplicities is None:
        multiplicities = {"adsorbate_slab": 1, "adsorbate": 1, "slab": 1}
    quantum_cluster = embedded_adsorbed_cluster[quantum_cluster_indices]
    ecp_region = embedded_adsorbed_cluster[ecp_region_indices]

    # Get the indices of the adsorbates from the quantum cluster
    adsorbate_indices = [
        i
        for i in range(len(quantum_cluster))
        if quantum_cluster.get_array("atom_type")[i] == "adsorbate"
    ]

    # Get the indices of the slab from the quantum cluster
    slab_indices = [
        i
        for i in range(len(quantum_cluster))
        if quantum_cluster.get_array("atom_type")[i] != "adsorbate"
    ]

    # Get the charge of the quantum cluster
    charge = int(sum(quantum_cluster.get_array("oxi_states")))

    # Get the total number of core electrons for the quantum cluster
    core = {
        "adsorbate_slab": sum(
            [element_info[atom.symbol]["core"] for atom in quantum_cluster]
        ),
        "adsorbate": sum(
            [
                element_info[atom.symbol]["core"]
                for atom in quantum_cluster
                if atom.index in adsorbate_indices
            ]
        ),
        "slab": sum(
            [
                element_info[atom.symbol]["core"]
                for atom in quantum_cluster
                if atom.index in slab_indices
            ]
        ),
    }

    # Create the coords strings for the adsorbate-slab complex, adsorbate, and slab
    coords_block = {
        "adsorbate_slab": f"""charge={charge}
mult={multiplicities['adsorbate_slab']}
core={int(core['adsorbate_slab']/2)}
unit=angs
geom=xyz
""",
        "adsorbate": f"""charge=0
mult={multiplicities['adsorbate']}
core={int(core['adsorbate']/2)}
unit=angs
geom=xyz
""",
        "slab": f"""charge={charge}
mult={multiplicities['slab']}
core={int(core['slab']/2)}
unit=angs
geom=xyz
""",
    }

    # Set the number of atoms for each system
    if include_cp:
        coords_block["adsorbate_slab"] += (
            f"{len(quantum_cluster) + len(ecp_region)}\n\n"
        )
        coords_block["adsorbate"] += f"{len(quantum_cluster)}\n\n"
        coords_block["slab"] += f"{len(quantum_cluster) + len(ecp_region)}\n\n"
    else:
        coords_block["adsorbate_slab"] += (
            f"{len(quantum_cluster) + len(ecp_region)}\n\n"
        )
        coords_block["adsorbate"] += f"{len(adsorbate_indices)}\n\n"
        coords_block["slab"] += f"{len(slab_indices) + len(ecp_region)}\n\n"

    for i, atom in enumerate(quantum_cluster):
        # Create the coords section for the adsorbate-slab complex
        coords_block["adsorbate_slab"] += create_atom_coord_string(atom=atom)

        # Create the coords section for the adsorbate and slab
        if i in adsorbate_indices:
            coords_block["adsorbate"] += create_atom_coord_string(atom=atom)
            if include_cp:
                coords_block["slab"] += create_atom_coord_string(atom=atom)
        elif i in slab_indices:
            coords_block["slab"] += create_atom_coord_string(atom=atom)
            if include_cp:
                coords_block["adsorbate"] += create_atom_coord_string(atom=atom)

    # Create the coords section for the ECP region
    for atom in ecp_region:
        coords_block["adsorbate_slab"] += create_atom_coord_string(atom=atom)
        coords_block["slab"] += create_atom_coord_string(atom=atom)

    # Adding the ghost atoms for the counterpoise correction
    for system in ["adsorbate_slab", "adsorbate", "slab"]:
        coords_block[system] += "\nghost=serialno\n"
        if include_cp and system in ["adsorbate"]:
            coords_block[system] += ",".join(
                [str(atom_idx + 1) for atom_idx in slab_indices]
            )
        elif include_cp and system in ["slab"]:
            coords_block[system] += ",".join(
                [str(atom_idx + 1) for atom_idx in adsorbate_indices]
            )
        coords_block[system] += "\n\n"

    return coords_block


def generate_mrcc_point_charge_block(
    embedded_adsorbed_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
) -> str:
    """
    Create the point charge block for the MRCC input file. This requires the embedded_cluster Atoms object containing both atom_type and oxi_states arrays, as well as the indices of the quantum cluster and ECP region. Such arrays are created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.

    Parameters
    ----------
    embedded_adsorbed_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file. This object can be created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    quantum_cluster_indices
        A list of lists containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.

    Returns
    -------
    str
        The point charge block for the MRCC input file.
    """

    # Get the oxi_states arrays from the embedded_cluster
    oxi_states = embedded_adsorbed_cluster.get_array("oxi_states")

    # Check that none of the indices in quantum_cluster_indices are in ecp_region_indices
    if not np.all([x not in ecp_region_indices for x in quantum_cluster_indices]):
        raise ValueError("An atom in the quantum cluster is also in the ECP region.")

    # Get the number of point charges for this system. There is a point charge associated with each capped ECP as well.
    pc_region_indices = ecp_region_indices + [
        atom.index
        for atom in embedded_adsorbed_cluster
        if atom not in quantum_cluster_indices + ecp_region_indices
    ]

    num_pc = len(pc_region_indices)
    pc_block = f"qmmm=Amber\npointcharges\n{num_pc}\n"

    # Add the ecp_region indices
    for i in pc_region_indices:
        position = embedded_adsorbed_cluster[i].position
        pc_block += f"  {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f} {oxi_states[i]:-16.11f}\n"

    return pc_block


def create_orca_eint_blocks(
    embedded_adsorbed_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
    element_info: dict[ElementStr, ElementInfo] | None = None,
    pal_nprocs_block: dict[str, int] | None = None,
    method_block: dict[str, str] | None = None,
    scf_block: dict[str, str] | None = None,
    ecp_info: dict[ElementStr, str] | None = None,
    include_cp: bool = True,
    multiplicities: MultiplicityDict | None = None,
) -> BlockInfo:
    """
    Creates the orcablocks input for the ORCA ASE calculator.

    Parameters
    ----------
    embedded_adsorbed_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This object is created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    quantum_cluster_indices
        A list containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list containing the indices of the atoms in each ECP region. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    pal_nprocs_block
        A dictionary with the number of processors for the PAL block as 'nprocs' and the maximum memory-per-core in megabytes blocks as 'maxcore'.
    method_block
        A dictionary that contains the method block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
    scf_block
        A dictionary that contains the SCF block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
    ecp_info
        A dictionary with the ECP data (in ORCA format) for the cations in the ECP region.
    include_cp
        If True, the coords strings will include the counterpoise correction for the adsorbate and slab.
    multiplicities
        The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.

    Returns
    -------
    BlockInfo
        The ORCA input block (to be put in 'orcablocks' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
    """

    # First generate the preamble block
    preamble_block = generate_orca_input_preamble(
        embedded_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=quantum_cluster_indices,
        element_info=element_info,
        pal_nprocs_block=pal_nprocs_block,
        method_block=method_block,
        scf_block=scf_block,
    )

    # Generate the coords block
    coords_block = generate_coords_block(
        embedded_adsorbed_cluster=embedded_adsorbed_cluster,
        quantum_cluster_indices=quantum_cluster_indices,
        ecp_region_indices=ecp_region_indices,
        ecp_info=ecp_info,
        include_cp=include_cp,
        multiplicities=multiplicities,
    )

    # Combine the blocks
    return {
        "adsorbate_slab": preamble_block + coords_block["adsorbate_slab"],
        "adsorbate": preamble_block + coords_block["adsorbate"],
        "slab": preamble_block + coords_block["slab"],
    }


def create_atom_coord_string(
    atom: Atom,
    is_ghost_atom: bool = False,
    atom_ecp_info: str | None = None,
    pc_charge: float | None = None,
) -> str:
    """
    Creates a string containing the Atom symbol and coordinates in the ORCA input file format, with additional information for atoms in the ECP region as well as ghost atoms.

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


def generate_coords_block(
    embedded_adsorbed_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
    ecp_info: dict[ElementStr, str] | None = None,
    include_cp: bool = True,
    multiplicities: MultiplicityDict | None = None,
) -> BlockInfo:
    """
    Generates the coordinates block for the ORCA input file. This includes the coordinates of the quantum cluster, the ECP region, and the point charges. It will return three strings for the adsorbate-slab complex, adsorbate and slab.

    Parameters
    ----------
    embedded_adsorbed_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file, as well as the atom type. This Atoms object is typically produced from [quacc.atoms.skzcam.create_skzcam_clusters][].
    quantum_cluster_indices
        A list containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list containing the indices of the atoms in each ECP region. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_info
        A dictionary with the ECP data (in ORCA format) for the cations in the ECP region.
    include_cp
        If True, the coords strings will include the counterpoise correction for the adsorbate and slab.
    multiplicities
        The multiplicity of the adsorbate-slab complex, adsorbate and slab respectively, with the keys 'adsorbate_slab', 'adsorbate', and 'slab'.

    Returns
    -------
    BlockInfo
        The coordinates block for the ORCA input block (to be put in 'orcablocks' parameter) as a string for the adsorbate-slab complex, the adsorbate, and the slab in a dictionary with the keys 'adsorbate_slab', 'adsorbate', and 'slab' respectively.
    """

    # Create the quantum cluster and ECP region cluster
    if multiplicities is None:
        multiplicities = {"adsorbate_slab": 1, "adsorbate": 1, "slab": 1}
    quantum_cluster = embedded_adsorbed_cluster[quantum_cluster_indices]
    ecp_region = embedded_adsorbed_cluster[ecp_region_indices]

    # Get the indices of the adsorbates from the quantum cluster
    adsorbate_indices = [
        i
        for i in range(len(quantum_cluster))
        if quantum_cluster.get_array("atom_type")[i] == "adsorbate"
    ]

    # Get the indices of the slab from the quantum cluster
    slab_indices = [
        i
        for i in range(len(quantum_cluster))
        if quantum_cluster.get_array("atom_type")[i] != "adsorbate"
    ]

    # Get the charge of the quantum cluster
    charge = int(sum(quantum_cluster.get_array("oxi_states")))

    # Create the coords strings for the adsorbate-slab complex, adsorbate, and slab
    coords_block = {
        "adsorbate_slab": f"""%coords
CTyp xyz
Mult {multiplicities['adsorbate_slab']}
Units angs
Charge {charge}
coords
""",
        "adsorbate": f"""%coords
CTyp xyz
Mult {multiplicities['adsorbate']}
Units angs
Charge 0
coords
""",
        "slab": f"""%coords
CTyp xyz
Mult {multiplicities['slab']}
Units angs
Charge {charge}
coords
""",
    }

    for i, atom in enumerate(quantum_cluster):
        # Create the coords section for the adsorbate-slab complex
        coords_block["adsorbate_slab"] += create_atom_coord_string(atom=atom)

        # Create the coords section for the adsorbate and slab
        if i in adsorbate_indices:
            coords_block["adsorbate"] += create_atom_coord_string(atom=atom)
            if include_cp:
                coords_block["slab"] += create_atom_coord_string(
                    atom=atom, is_ghost_atom=True
                )
        elif i in slab_indices:
            coords_block["slab"] += create_atom_coord_string(atom=atom)
            if include_cp:
                coords_block["adsorbate"] += create_atom_coord_string(
                    atom=atom, is_ghost_atom=True
                )

    # Create the coords section for the ECP region
    ecp_region_coords_section = ""
    for i, atom in enumerate(ecp_region):
        atom_ecp_info = format_ecp_info(atom_ecp_info=ecp_info[atom.symbol])
        ecp_region_coords_section += create_atom_coord_string(
            atom=atom,
            atom_ecp_info=atom_ecp_info,
            pc_charge=ecp_region.get_array("oxi_states")[i],
        )

    # Add the ECP region coords section to the ads_slab_coords string
    coords_block["adsorbate_slab"] += f"{ecp_region_coords_section}end\nend\n"
    coords_block["slab"] += f"{ecp_region_coords_section}end\nend\n"
    coords_block["adsorbate"] += "end\nend\n"

    return coords_block


def format_ecp_info(atom_ecp_info: str) -> str:
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


def generate_orca_input_preamble(
    embedded_cluster: Atoms,
    quantum_cluster_indices: list[int],
    element_info: dict[ElementStr, ElementInfo] | None = None,
    pal_nprocs_block: dict[str, int] | None = None,
    method_block: dict[str, str] | None = None,
    scf_block: dict[str, str] | None = None,
) -> str:
    """
    From the quantum cluster Atoms object, generate the ORCA input preamble for the basis, method, pal, and scf blocks.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun ChemShell file. This object is created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    quantum_cluster_indices
        A list containing the indices of the atoms of embedded_cluster that form a quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    element_info
        A dictionary with elements as keys which gives the (1) number of core electrons as 'core', (2) basis set as 'basis', (3) effective core potential as 'ecp', (4) resolution-of-identity/density-fitting auxiliary basis set for DFT/HF calculations as 'ri_scf_basis' and (5) resolution-of-identity/density-fitting for correlated wave-function methods as 'ri_cwft_basis'.
    pal_nprocs_block
        A dictionary with the number of processors for the PAL block as 'nprocs' and the maximum memory-per-core blocks as 'maxcore'.
    method_block
        A dictionary that contains the method block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.
    scf_block
        A dictionary that contains the SCF block for the ORCA input file. The key is the ORCA setting and the value is that setting's value.

    Returns
    -------
    str
        The ORCA input preamble.
    """

    # Create the quantum cluster
    quantum_cluster = embedded_cluster[quantum_cluster_indices]

    # Get the set of element symbols from the quantum cluster
    element_symbols = list(set(quantum_cluster.get_chemical_symbols()))
    element_symbols.sort()

    # Check all element symbols are provided in element_info keys
    if element_info is not None and not all(
        element in element_info for element in element_symbols
    ):
        raise ValueError(
            "Not all element symbols are provided in the element_info dictionary."
        )

    # Initialize preamble_info
    preamble_info = """"""

    # Add the pal_nprocs_block
    if pal_nprocs_block is not None:
        preamble_info += f"%pal nprocs {pal_nprocs_block['nprocs']} end\n"
        preamble_info += f"%maxcore {pal_nprocs_block['maxcore']} end\n"

    # Add pointcharge file to read. It will be assumed that it is in the same folder as the input file
    preamble_info += '%pointcharges "orca.pc"\n'

    # Make the method block
    if method_block is not None and element_info is not None:
        preamble_info += "%method\n"
    # Iterate through the keys of method_block and add key value
    if method_block is not None:
        for key in method_block:
            preamble_info += f"{key} {method_block[key]}\n"
    # Iterate over the core value for each element (if it has been given)
    if element_info is not None:
        for element in element_symbols:
            if "core" in element_info[element]:
                preamble_info += (
                    f"NewNCore {element} {element_info[element]['core']} end\n"
                )
    if method_block is not None and element_info is not None:
        preamble_info += "end\n"

    # Make the basis block

    # First check if the basis key is the same for all elements. We use """ here because an option for these keys is "AutoAux"
    if element_info is not None:
        preamble_info += "%basis\n"
        if len({element_info[element]["basis"] for element in element_symbols}) == 1:
            preamble_info += f"""Basis {element_info[element_symbols[0]]['basis']}\n"""
        else:
            for element in element_symbols:
                element_basis = element_info[element]["basis"]
                preamble_info += f"""NewGTO {element} "{element_basis}" end\n"""

        # Do the same for ri_scf_basis and ri_cwft_basis.
        if (
            len({element_info[element]["ri_scf_basis"] for element in element_symbols})
            == 1
        ):
            preamble_info += (
                f"""Aux {element_info[element_symbols[0]]['ri_scf_basis']}\n"""
            )
        else:
            for element in element_symbols:
                element_basis = element_info[element]["ri_scf_basis"]
                preamble_info += f'NewAuxJGTO {element} "{element_basis}" end\n'

        if (
            len(
                list(
                    {
                        element_info[element]["ri_cwft_basis"]
                        for element in element_symbols
                    }
                )
            )
            == 1
        ):
            preamble_info += (
                f"""AuxC {element_info[element_symbols[0]]['ri_cwft_basis']}\n"""
            )
        else:
            for element in element_symbols:
                element_basis = element_info[element]["ri_cwft_basis"]
                preamble_info += f"""NewAuxCGTO {element} "{element_basis}" end\n"""

        preamble_info += "end\n"

    # Write the scf block
    if scf_block is not None:
        preamble_info += "%scf\n"
        for key in scf_block:
            preamble_info += f"""{key} {scf_block[key]}\n"""
        preamble_info += "end\n"

    return preamble_info


def create_orca_point_charge_file(
    embedded_cluster: Atoms,
    quantum_cluster_indices: list[int],
    ecp_region_indices: list[int],
    pc_file: str | Path,
) -> None:
    """
    Create a point charge file that can be read by ORCA. This requires the embedded_cluster Atoms object containing both atom_type and oxi_states arrays, as well as the indices of the quantum cluster and ECP region.

    Parameters
    ----------
    embedded_cluster
        The ASE Atoms object containing the atomic coordinates and atomic charges from the .pun file. This object is created by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    quantum_cluster_indices
        A list of lists containing the indices of the atoms in each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    ecp_region_indices
        A list of lists containing the indices of the atoms in the ECP region for each quantum cluster. These indices are provided by the [quacc.atoms.skzcam.create_skzcam_clusters][] function.
    pc_file
        A file containing the point charges to be written by ORCA.

    Returns
    -------
    None
    """

    # Get the oxi_states arrays from the embedded_cluster
    oxi_states = embedded_cluster.get_array("oxi_states")

    # Check that none of the indices in quantum_cluster_indices are in ecp_region_indices
    if not np.all([x not in ecp_region_indices for x in quantum_cluster_indices]):
        raise ValueError("An atom in the quantum cluster is also in the ECP region.")

    # Get the number of point charges for this system
    total_indices = quantum_cluster_indices + ecp_region_indices
    num_pc = len(embedded_cluster) - len(total_indices)
    counter = 0
    with Path.open(pc_file, "w") as f:
        # Write the number of point charges first
        f.write(f"{num_pc}\n")
        for i in range(len(embedded_cluster)):
            if i not in total_indices:
                counter += 1
                position = embedded_cluster[i].position
                if counter != num_pc:
                    f.write(
                        f"{oxi_states[i]:-16.11f} {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f}\n"
                    )
                else:
                    f.write(
                        f"{oxi_states[i]:-16.11f} {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f}"
                    )


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
        The indices of the atoms that are at the 'center' of the slab right beneath the adsorbate.

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
    slab_center_idx = next(
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
        slab_center_idx,
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
    chemsh_embedded_cluster.save(filename=Path(filepath).with_suffix(".pun"), fmt="pun")

    if write_xyz_file:
        # XYZ for visualisation
        chemsh_embedded_cluster.save(
            filename=Path(filepath).with_suffix(".xyz"), fmt="xyz"
        )


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
    embedded_cluster = convert_pun_to_atoms(
        pun_file=pun_file, atom_oxi_states=atom_oxi_states
    )

    # Get distances of all atoms from the cluster center
    atom_center_distances = _get_atom_distances(
        embedded_cluster=embedded_cluster, center_position=center_position
    )

    # Determine the cation shells from the center of the embedded cluster
    _, cation_shells_idx = _find_cation_shells(
        embedded_cluster=embedded_cluster,
        distances=atom_center_distances,
        shell_width=shell_width,
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
        embedded_cluster=embedded_cluster,
        quantum_cluster_indices=quantum_cluster_indices,
        dist_matrix=embedded_cluster_all_dist,
        ecp_dist=ecp_dist,
    )

    # Write the quantum clusters to files
    if write_clusters:
        for idx in range(len(quantum_cluster_indices)):
            quantum_atoms = embedded_cluster[quantum_cluster_indices[idx]]
            if write_include_ecp:
                ecp_atoms = embedded_cluster[ecp_region_indices[idx]]
                ecp_atoms.set_chemical_symbols(np.array(["U"] * len(ecp_atoms)))
                cluster_atoms = quantum_atoms + ecp_atoms
            else:
                cluster_atoms = quantum_atoms
            write(Path(write_clusters_path, f"SKZCAM_cluster_{idx}.xyz"), cluster_atoms)

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
    ecp_region_indices: list[list[int]] | None = None,
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
    ecp_region_indices
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
            list(range(len(adsorbate))) + [idx + len(adsorbate) for idx in cluster]
            for cluster in quantum_cluster_indices
        ]
    if ecp_region_indices is not None:
        ecp_region_indices = [
            [idx + len(adsorbate) for idx in cluster] for cluster in ecp_region_indices
        ]

    return embedded_adsorbate_cluster, quantum_cluster_indices, ecp_region_indices


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
    embedded_cluster: Atoms, distances: NDArray, shell_width: float = 0.1
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

    current_point = distances_sorted[0]
    current_shell = [current_point]
    current_shell_idx = [distances_sorted_indices[0]]

    for idx, point in enumerate(distances_sorted[1:]):
        if point <= current_point + shell_width:
            current_shell.append(point)
            current_shell_idx.append(distances_sorted_indices[idx + 1])
        else:
            shells.append(current_shell)
            shells_indices.append(current_shell_idx)
            current_shell = [point]
            current_shell_idx = [distances_sorted_indices[idx + 1]]
        current_point = point
    shells.append(current_shell)
    shells_indices.append(current_shell_idx)

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
