from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

from ase.io import read
from ase.units import Hartree

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import MRCCEnergyInfo


def read_geom_mrccinp(file_path: Path | str) -> Atoms:
    """
    Read geometry from an MRCC input file.

    Parameters
    ----------
    file_path: Path | str
        The path to the MRCC input file.

    Returns
    -------
    Atoms
        Atoms object with the geometry.
    """

    # Get the lines as a list
    with Path.open(file_path) as fd:
        lines = fd.readlines()
    xyz_line_index = [index for index, line in enumerate(lines) if "xyz" in line]

    if len(xyz_line_index) != 1:
        raise ValueError("Geometry incorrectly provided in MRCC input file")

    xyz_line_index = xyz_line_index[0]

    # Get the number of atoms
    atoms_length = int(lines[xyz_line_index + 1])

    # Format and send the string to be read by ase.io.read()
    xyz_text = f"{atoms_length}\n geometry\n"
    for line in lines[xyz_line_index + 3 : xyz_line_index + 3 + atoms_length]:
        xyz_text += line
    atoms = read(StringIO(xyz_text), format="xyz")

    # Remove PBC and set the unit cell to zero as MRCC is a molecular code.
    atoms.pbc = False
    atoms.set_cell([0.0, 0.0, 0.0])

    return atoms


def write_mrcc(file_path: Path | str, atoms: Atoms, parameters: dict[str, str]) -> None:
    """
    Write MRCC input file given the Atoms object and the parameters.

    Parameters
    ----------
    file_path : Path | str
        File path to write the MRCC input file.
    atoms : Atoms
        Atoms object with the geometry.
    parameters : dict[str,str]
        Dictionary with the parameters to be written in the MRCC input file. The keys are the input keyword and the values are the input values.
    """

    with Path.open(file_path, "w") as file_path:
        # Write the MRCC input file
        for key, value in parameters.items():
            file_path.write(f"{key}={value}\n")

        if "geom" not in parameters:
            # If the geometry is not provided in the MRCC blocks, write it here.
            ghost_list = []  # List of indices of the ghost atoms.
            file_path.write(f"geom=xyz\n{len(atoms)}\n\n")
            for atom_idx, atom in enumerate(atoms):
                if atom.tag == 71:  # 71 is ascii G (Ghost)
                    ghost_list += [atom_idx + 1]

                symbol = atom.symbol
                position = atom.position
                file_path.write(
                    f"{symbol.ljust(3)} {position[0]:-16.11f} {position[1]:-16.11f} {position[2]:-16.11f}\n"
                )

            if ghost_list and "ghost" not in parameters:
                file_path.write("\nghost=serialno\n")
                file_path.write(",".join([str(atom_idx) for atom_idx in ghost_list]))


def read_energy(lines: list[str]) -> MRCCEnergyInfo:
    """
    Reads the energy components (SCF energy, MP2 correlation energy, CCSD correlation energy, CCSD(T) correlation energy) from the MRCC output file where available.

    Parameters
    ----------
    lines : list[str]
        List of lines read from the MRCC output file.

    Returns
    -------
    MRCCEnergyInfo
        Dictionary with the energy components. The keys are the following:
        - energy : float <-- Total energy which will not be computed in this function.
        - scf_energy : float <-- SCF energy.
        - mp2_corr_energy : float <-- MP2 correlation energy.
        - ccsd_corr_energy : float <-- CCSD correlation energy.
        - ccsdt_corr_energy : float <-- CCSD(T) correlation energy.
    """

    energy_dict = {
        "energy": None,
        "scf_energy": None,
        "mp2_corr_energy": None,
        "ccsd_corr_energy": None,
        "ccsdt_corr_energy": None,
    }

    for line in lines:
        if "FINAL HARTREE-FOCK ENERGY" in line or "FINAL KOHN-SHAM ENERGY" in line:
            energy_dict["scf_energy"] = float(line.split()[-2]) * Hartree
        elif "MP2 correlation energy" in line:
            energy_dict["mp2_corr_energy"] = float(line.split()[-1]) * Hartree
        elif "CCSD correlation energy" in line:
            energy_dict["ccsd_corr_energy"] = float(line.split()[-1]) * Hartree
        elif "CCSD(T) correlation energy" in line:
            energy_dict["ccsdt_corr_energy"] = float(line.split()[-1]) * Hartree

    return energy_dict


def read_mrcc_outputs(output_file_path: Path | str) -> MRCCEnergyInfo:
    """
    Reads the energy components (SCF energy, MP2 correlation energy, CCSD correlation energy, CCSD(T) correlation energy) from the MRCC output file where available and calculates the total energy (based on the highest level of theory)

    Parameters
    ----------
    output_file_path : Path | str
        Path to the MRCC output file.

    Returns
    -------
    MRCCEnergyInfo
        Dictionary with the energy components. The keys are the following:
        - energy : float | None <-- Total energy of highest available level.
        - scf_energy : float | None <-- SCF energy.
        - mp2_corr_energy : float | None <-- MP2 correlation energy.
        - ccsd_corr_energy : float | None <-- CCSD correlation energy.
        - ccsdt_corr_energy : float | None <-- CCSD(T) correlation energy.
    """
    with Path.open(output_file_path) as output_textio:
        lines = output_textio.readlines()

    energy_dict = read_energy(lines)

    # Raise error if scf_energy is None
    if energy_dict["scf_energy"] is None:
        raise ValueError("SCF energy not found in MRCC output file")

    if energy_dict["ccsdt_corr_energy"] is not None:
        energy_dict["energy"] = (
            energy_dict["scf_energy"] + energy_dict["ccsdt_corr_energy"]
        )
    elif energy_dict["ccsd_corr_energy"] is not None:
        energy_dict["energy"] = (
            energy_dict["scf_energy"] + energy_dict["ccsd_corr_energy"]
        )
    elif energy_dict["mp2_corr_energy"] is not None:
        energy_dict["energy"] = (
            energy_dict["scf_energy"] + energy_dict["mp2_corr_energy"]
        )
    else:
        energy_dict["energy"] = energy_dict["scf_energy"]

    return energy_dict
