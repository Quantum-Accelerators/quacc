import re
from io import StringIO
from pathlib import Path
from typing import List, Optional

import numpy as np

from ase.io import read
from ase.units import Bohr, Hartree
from ase.utils import reader, writer

# Made from NWChem interface

@reader
def read_geom_mrccinp(fd):
    """Method to read geometry from an MRCC input file."""
    lines = fd.readlines()
    atomslength = 0
    atomslenline = -1

    # Find geometry region of input file.
    for index, line in enumerate(lines):
        if 'xyz' in line:
            startline = index + 3
            atomslenline = index + 1
        elif index == atomslenline:
            atomslength = int(line)
            break

    # Format and send to read_xyz.
    xyz_text = '%i\n' % (atomslength)
    xyz_text += ' geometry\n'
    for line in lines[startline:startline+atomslength]:
        xyz_text += line
    atoms = read(StringIO(xyz_text), format='xyz')
    atoms.set_cell((0., 0., 0.))  # no unit cell defined

    return atoms


@writer
def write_mrcc(fd, atoms, params):
    # List of indices of the ghost atoms.
    ghost_list = []

    for key, value in params['mrccinput'].items():
        fd.write(f"{key}={value}\n")
    
    fd.write(f"{params['mrccblocks']} \n")

    if 'geom' not in params['mrccblocks']:
        fd.write("charge=%d\n" % params['charge'])
        fd.write("mult=%d\n" % params['mult'])
        fd.write('geom=xyz\n')
        fd.write("%d\n\n" % len(atoms))
        for atom_idx, atom in enumerate(atoms):
            if atom.tag == 71:  # 71 is ascii G (Ghost)
                ghost_list += [atom_idx+1]
                symbol = atom.symbol + '   '
            else:
                symbol = atom.symbol + '   '
            fd.write(
                symbol
                + str(atom.position[0])
                + " "
                + str(atom.position[1])
                + " "
                + str(atom.position[2])
                + "\n"
            )

        if len(ghost_list) > 0:
            fd.write("\nghost=serialno")
            fd.write(','.join([str(atom_idx) for atom_idx in ghost_list]))


def read_energy(lines: List[str]) -> Optional[float]:
    """Read energy."""
    scf_energy = None
    mp2_corr_energy = None
    ccsd_corr_energy = None
    ccsdt_corr_energy = None
    for line in lines:
        if 'FINAL HARTREE-FOCK ENERGY' in line:
            scf_energy = float(line.split()[-2])*Hartree
        elif 'FINAL KOHN-SHAM ENERGY' in line:
            scf_energy = float(line.split()[-2])*Hartree
        elif 'MP2 correlation energy' in line:
            mp2_corr_energy = float(line.split()[-1])*Hartree
        elif 'CCSD correlation energy' in line:
            ccsd_corr_energy = float(line.split()[-1])*Hartree
        elif 'CCSD(T) correlation energy' in line:
            ccsdt_corr_energy = float(line.split()[-1])*Hartree

    return scf_energy, mp2_corr_energy, ccsd_corr_energy, ccsdt_corr_energy

@reader
def read_mrcc_outputs(fd):
    """ From the MRCC output file: Read Energy
    """
    lines = fd.readlines()

    scf_energy, mp2_corr_energy, ccsd_corr_energy, ccsdt_corr_energy = read_energy(lines)

    results = {}

    if scf_energy is not None:
        results['scf_energy'] = scf_energy
        results['energy'] = scf_energy
    if mp2_corr_energy is not None:
        results['mp2_corr_energy'] = mp2_corr_energy
        results['energy'] = scf_energy + mp2_corr_energy
    if ccsd_corr_energy is not None:
        results['ccsd_corr_energy'] = ccsd_corr_energy
        results['energy'] = scf_energy + ccsd_corr_energy
    if ccsdt_corr_energy is not None:
        results['ccsd(t)_corr_energy'] = ccsdt_corr_energy
        results['energy'] = scf_energy + ccsdt_corr_energy

    return results




