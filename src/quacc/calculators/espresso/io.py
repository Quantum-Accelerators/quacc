"""I/O utilities for the espresso calculator."""
# 1. use the globals() dictionary
# 2. use sys.modules[__name__] and getattr
# 3. keep track of all functions in a dictionary
from __future__ import annotations

import re

import numpy as np
from ase import Atoms
from ase.io import read as _read
from ase.io import write as _write
from ase.units import Bohr

from quacc.calculators.espresso.utils import namelist_to_string

# Three ways to call a function based on a string:
freg = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?")


def write(filename, atoms, format="pw", properties=None, **kwargs):
    with open(filename, "w") as fd:
        write_espresso_dict[format](fd, atoms=atoms, properties=properties, **kwargs)
    # lines = write_namelist(filename, parameters)
    # specific_lines = func_dict[format][1](**kwargs)
    # lines += specific_lines
    # with open(filename, 'w') as fd:
    #    fd.write(''.join(lines))


def read(filename, format="pw"):
    with open(filename) as fd:
        return read_espresso_dict[format](fd)


def write_espresso_io(fd, **kwargs):
    "For simple binaries, this is enough."
    pwi = namelist_to_string(kwargs["input_data"])
    fd.write("".join(pwi))


def write_espresso_ph(fd, **kwargs):
    input_data = kwargs.get("input_data", None)
    qpts = kwargs.get("qpts", None)
    pwi = namelist_to_string(input_data)[:-1]
    fd.write("".join(pwi))
    qplot = input_data["inputph"].get("qplot", False)
    ldisp = input_data["inputph"].get("ldisp", False)
    if qplot:
        fd.write(f"{len(qpts)}\n")
        for qpt in qpts:
            fd.write("{0:8.4f} {1:8.4f} {2:8.4f}\n".format(*qpt))
    elif not ldisp:
        fd.write("0.0 0.0 0.0\n")


def read_espresso_ph(fd):
    QPOINTS = r"(?i)^\s*Calculation\s*of\s*q"
    NKPTS = r"(?i)^\s*number\s*of\s*k\s*points\s*"
    DIEL = r"(?i)^\s*Dielectric\s*constant\s*in\s*cartesian\s*axis\s*$"
    BORN = r"(?i)^\s*Effective\s*charges\s*\(d\s*Force\s*/\s*dE\)"
    POLA = r"(?i)^\s*Polarizability\s*(a.u.)\^3"
    MODE = r"(?i)^\s*(mode\s*#\s*\d\s*)+"
    EQPOINTS = r"(?i)^\s*Number\s*of\s*q\s*in\s*the\s*star\s*=\s*"
    DIAG = r"(?i)^\s*Diagonalizing\s*the\s*dynamical\s*matrix\s*$"
    MODE_SYM = r"(?i)^\s*Mode\s*symmetry,\s*"
    BORN_DFPT = r"(?i)^\s*Effective\s*charges\s*\(d\s*P\s*/\s*du\)"
    POSITIONS = r"(?i)^\s*site\s*n\..*\(alat\s*units\)"
    ALAT = r"(?i)^\s*celldm\(1\)="
    CELL = r"^\s*crystal\s*axes:\s*\(cart.\s*coord.\s*in\s*units\s*of\s*alat\)"

    output = {
        QPOINTS: [],
        NKPTS: [],
        DIEL: [],
        BORN: [],
        BORN_DFPT: [],
        POLA: [],
        MODE: [],
        EQPOINTS: [],
        DIAG: [],
        MODE_SYM: [],
        POSITIONS: [],
        ALAT: [],
        CELL: [],
    }

    names = {
        QPOINTS: "qpoints",
        NKPTS: "kpoints",
        DIEL: "dieltensor",
        BORN: "borneffcharge",
        BORN_DFPT: "borneffcharge_dfpt",
        POLA: "polarizability",
        MODE: "modes",
        EQPOINTS: "eqpoints",
        DIAG: "freqs",
        MODE_SYM: "mode_symmetries",
        POSITIONS: "positions",
        ALAT: "alat",
        CELL: "cell",
    }

    results = {}
    fdo_lines = [i for i in fd.read().splitlines() if i]
    n_lines = len(fdo_lines)

    for idx, line in enumerate(fdo_lines):
        for key in output:
            if bool(re.match(key, line)):
                output[key].append(idx)

    output = {key: np.array(value) for key, value in output.items()}

    def _read_qpoints(idx):
        match = re.findall(freg, fdo_lines[idx])
        return tuple([float(x) for x in match])

    def _read_kpoints(idx):
        n_kpts = int(re.findall(freg, fdo_lines[idx])[0])
        kpts = []
        for line in fdo_lines[idx + 2 : idx + 2 + n_kpts]:
            if bool(re.search(r"^\s*k\(.*wk", line)):
                kpts.append([float(x) for x in re.findall(freg, line)[1:]])
        return np.array(kpts)

    def _read_modes(idx):
        n = 1
        n_modes = len(re.findall(r"mode", fdo_lines[idx]))
        modes = []
        while len(modes) == 0 or bool(re.match(r"^\s*\(", fdo_lines[idx + n])):
            tmp = re.findall(freg, fdo_lines[idx + n])
            modes.append([float(x) for x in tmp])
            n += 1
        return np.hsplit(np.array(modes), n_modes)

    def _read_eqpoints(idx):
        n_star = int(re.findall(freg, fdo_lines[idx])[0])
        eqpoints = np.loadtxt(fdo_lines[idx + 2 : idx + 2 + n_star], usecols=(1, 2, 3))
        return eqpoints

    def _read_freqs(idx):
        n = 0
        freqs = []
        stop = 0
        while len(freqs) == 0 or stop < 2:
            if bool(re.search(r"^\s*freq", fdo_lines[idx + n])):
                tmp = re.findall(freg, fdo_lines[idx + n])[1]
                freqs.append(float(tmp))
            if bool(re.search(r"\*{5,}", fdo_lines[idx + n])):
                stop += 1
            n += 1
        return np.array(freqs)

    def _read_sym(idx):
        n = 1
        sym = {}
        while bool(re.match(r"^\s*freq", fdo_lines[idx + n])):
            r = re.findall("\\d+", fdo_lines[idx + n])
            r = tuple(range(int(r[0]), int(r[1]) + 1))
            sym[r] = fdo_lines[idx + n].split("-->")[1].strip()
            sym[r] = re.sub(r"\s+", " ", sym[r])
            n += 1
        return sym

    def _read_epsil(idx):
        epsil = np.zeros((3, 3))
        for n in range(1, 4):
            tmp = re.findall(freg, fdo_lines[idx + n])
            epsil[n - 1] = [float(x) for x in tmp]
        return epsil

    def _read_born(idx):
        n = 1
        born = []
        while idx + n < n_lines:
            if re.search(r"^\s*atom\s*\d\s*\S", fdo_lines[idx + n]):
                pass
            elif re.search(r"^\s*E(x|y|z)\s*\(", fdo_lines[idx + n]):
                tmp = re.findall(freg, fdo_lines[idx + n])
                born.append([float(x) for x in tmp])
            else:
                break
            n += 1
        born = np.array(born)
        return np.vsplit(born, len(born) // 3)

    def _read_born_dfpt(idx):
        n = 1
        born = []
        while idx + n < n_lines:
            if re.search(r"^\s*atom\s*\d\s*\S", fdo_lines[idx + n]):
                pass
            elif re.search(r"^\s*P(x|y|z)\s*\(", fdo_lines[idx + n]):
                tmp = re.findall(freg, fdo_lines[idx + n])
                born.append([float(x) for x in tmp])
            else:
                break
            n += 1
        born = np.array(born)
        return np.vsplit(born, len(born) // 3)

    def _read_pola(idx):
        pola = np.zeros((3, 3))
        for n in range(1, 4):
            tmp = re.findall(freg, fdo_lines[idx + n])[:3]
            pola[n - 1] = [float(x) for x in tmp]
        return pola

    def _read_positions(idx):
        positions = []
        symbols = []
        n = 1
        while re.findall(r"^\s*\d+", fdo_lines[idx + n]):
            symbols.append(fdo_lines[idx + n].split()[1])
            positions.append(
                [float(x) for x in re.findall(freg, fdo_lines[idx + n])[-3:]]
            )
            n += 1
        atoms = Atoms(positions=positions, symbols=symbols)
        atoms.pbc = True
        return atoms

    def _read_alat(idx):
        return float(re.findall(freg, fdo_lines[idx])[1])

    def _read_cell(idx):
        cell = []
        n = 1
        while re.findall(r"^\s*a\(\d\)", fdo_lines[idx + n]):
            cell.append([float(x) for x in re.findall(freg, fdo_lines[idx + n])[-3:]])
            n += 1
        return np.array(cell)

    properties = {
        NKPTS: _read_kpoints,
        DIEL: _read_epsil,
        BORN: _read_born,
        BORN_DFPT: _read_born_dfpt,
        POLA: _read_pola,
        MODE: _read_modes,
        EQPOINTS: _read_eqpoints,
        DIAG: _read_freqs,
        MODE_SYM: _read_sym,
        POSITIONS: _read_positions,
        ALAT: _read_alat,
        CELL: _read_cell,
    }

    # Should this be printed only once?
    unique = {
        QPOINTS: True,
        NKPTS: False,
        DIEL: True,
        BORN: True,
        BORN_DFPT: True,
        POLA: True,
        MODE: False,
        EQPOINTS: True,
        DIAG: True,
        MODE_SYM: True,
        POSITIONS: True,
        ALAT: True,
        CELL: True,
    }

    # functions is a dict of functions to be applied to the output

    iblocks = np.append(output[QPOINTS], n_lines)

    for past, future in zip(iblocks[:-1], iblocks[1:]):
        qpoint = _read_qpoints(past)
        results[qpoint] = {}
        for prop in properties:
            p = (past < output[prop]) & (output[prop] < future)
            selected = output[prop][p]
            if len(selected) == 0:
                continue
            elif unique[prop]:
                idx = output[prop][p][-1]
                results[qpoint][names[prop]] = properties[prop](idx)
            else:
                tmp = {k + 1: 0 for k in range(len(selected))}
                for k, idx in enumerate(selected):
                    tmp[k + 1] = properties[prop](idx)
                results[qpoint][names[prop]] = tmp
        alat = results[qpoint].pop("alat", 1.0)
        atoms = results[qpoint].pop("positions", None)
        cell = results[qpoint].pop("cell", np.eye(3))
        if atoms:
            atoms.positions *= alat * Bohr
            atoms.cell = cell * alat * Bohr
            atoms.wrap()
            results[qpoint]["atoms"] = atoms
        results[qpoint]["atoms"]

    return results


def read_espresso_pw(filename):
    return _read(filename, format="espresso-out")


def write_espresso_pw(filename, atoms, properties, **kwargs):
    _write(filename, atoms, format="espresso-in", properties=properties, **kwargs)


read_espresso_dict = {"ph": read_espresso_ph, "pw": read_espresso_pw}

write_espresso_dict = {"ph": write_espresso_ph, "pw": write_espresso_pw}
