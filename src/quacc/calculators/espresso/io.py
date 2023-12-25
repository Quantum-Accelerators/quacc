"""I/O utilities for the espresso calculator."""
from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import Atoms
from ase.io import read as _read
from ase.io import write as _write
from ase.io.espresso import namelist_to_string
from ase.units import Bohr

if TYPE_CHECKING:
    from typing import Any, TextIO

freg = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?")


def write(
    filename: str | Path,
    atoms: Atoms,
    binary: str = "pw",
    properties: dict[str, Any] | None = None,
    **kwargs,
) -> None:
    """
    Wrapper around ase.io.espresso.write to read additional binary, such as ph.x.
    This is probably temporary until ase.io.espresso is updated.
    """

    filename_path = Path(filename).expanduser()
    with Path.open(filename_path, "w") as fd:
        write_espresso_dict[binary](fd, atoms=atoms, properties=properties, **kwargs)


def read(filename: str | Path, binary: str = "pw") -> dict[str, Any]:
    """
    Wrapper around ase.io.espresso.read to read additional binary, such as ph.x.
    This is probably temporary until ase.io.espresso is updated.
    """

    filename_path = Path(filename).expanduser()
    with Path.open(filename_path) as fd:
        return read_espresso_dict[binary](fd)


def write_espresso_io(fd: TextIO, **kwargs) -> None:
    """
    Function which writes input for simple espresso binaries without special cards.
    Non-exhaustive list (to complete):

    - pp.x
    - dynmat.x
    - projwfc.x

    Note: "EOF" is appended at the end of the file.
    (https://lists.quantum-espresso.org/pipermail/users/2020-November/046269.html)

    Parameters
    ----------
    fd
        The file descriptor of the input file.
    kwargs
        kwargs dictionary possibly containing the following keys:

        - input_data: dict

    Returns
    -------
    None
    """

    pwi = namelist_to_string(kwargs["input_data"])[:-1]
    pwi.append("EOF")
    fd.write("".join(pwi))


def write_espresso_ph(fd: TextIO, **kwargs) -> None:
    """
    Function that write the input file for a ph.x calculation. Normal namelist
    cards are passed in the input_data dictionary. Which can be either nested
    or flat, ASE style. The q-points are passed in the qpts list. If qplot is
    set to True then qpts is expected to be a list of list|tuple of length 4.
    Where the first three elements are the coordinates of the q-point in units
    of 2pi/alat and the last element is the weight of the q-point. if qplot is
    set to False then qpts is expected to be a simple list of length 4 (single
    q-point). Finally if ldisp is set to True, the above is discarded and the
    q-points are read from the nq1, nq2, nq3 cards in the input_data dictionary.

    Additionally, a nat_todo kwargs (list[int]) can be specified in the kwargs.
    It will be used if nat_todo is set to True in the input_data dictionary.

    Globally, this function follows the convention set in the ph.x documentation
    (https://www.quantum-espresso.org/Doc/INPUT_PH.html)

    Parameters
    ----------
    fd
        The file descriptor of the input file.
    kwargs
        kwargs dictionary possibly containing the following keys:

        - input_data: dict
        - qpts: list[list[float]] | list[tuple[float]] | list[float]
        - nat_todo: list[int]

    Returns
    -------
    None
    """

    input_data = kwargs.get("input_data")
    pwi = namelist_to_string(input_data)[:-1]
    fd.write("".join(pwi))

    qpts = kwargs.get("qpts")
    qplot = input_data["inputph"].get("qplot", False)
    ldisp = input_data["inputph"].get("ldisp", False)
    nat_todo = input_data["inputph"].get("nat_todo", 0)

    if qplot and ldisp:
        fd.write(f"{len(qpts)}\n")
        for qpt in qpts:
            fd.write(f"{qpt[0]:0.8f} {qpt[1]:0.8f} {qpt[2]:0.8f} {qpt[3]:1d}\n")
    else:
        fd.write(f"{qpts[0]:0.8f} {qpts[1]:0.8f} {qpts[2]:0.8f}\n")
    if nat_todo:
        nat_todo_index = kwargs.get("nat_todo")
        tmp = [str(i) for i in nat_todo_index]
        fd.write(" ".join(tmp))
        fd.write("\n")


def read_espresso_ph(fd: TextIO) -> dict[str, Any]:
    """
    Function that reads the output of a ph.x calculation. It returns a dictionary
    where each q-point is a key and the value is a dictionary with the following
    keys if available:

    - qpoints: The q-point in cartesian coordinates.
    - kpoints: The k-points in cartesian coordinates.
    - dieltensor: The dielectric tensor.
    - borneffcharge: The effective Born charges.
    - borneffcharge_dfpt: The effective Born charges from DFPT.
    - polarizability: The polarizability tensor.
    - modes: The phonon modes.
    - eqpoints: The symmetrically equivalent q-points.
    - freqs: The phonon frequencies.
    - mode_symmetries: The symmetries of the modes.
    - atoms: The atoms object.

    Some notes:

        - For some reason, the cell is not defined to high level of precision with
            ph.x. Be careful when using the atoms object retrieved from this function.
        - This function can be called on incomplete calculations i.e. if the calculation
            couldn't diagonalize the dynamical matrix for some q-points, the results for
            the other q-points will still be returned.

    Parameters
    ----------
    fd
        The file descriptor of the output file.

    Returns
    -------
    dict
        The results dictionnary as described above.
    """
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
    ELECTRON_PHONON = r"(?i)^\s*electron-phonon\s*interaction\s*...\s*$"

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
        ELECTRON_PHONON: [],
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
        ELECTRON_PHONON: "ep_data",
    }

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
        ELECTRON_PHONON: True,
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
        return tuple(float(x) for x in match)

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
        while not modes or bool(re.match(r"^\s*\(", fdo_lines[idx + n])):
            tmp = re.findall(freg, fdo_lines[idx + n])
            modes.append([float(x) for x in tmp])
            n += 1
        return np.hsplit(np.array(modes), n_modes)

    def _read_eqpoints(idx):
        n_star = int(re.findall(freg, fdo_lines[idx])[0])
        return np.loadtxt(
            fdo_lines[idx + 2 : idx + 2 + n_star], usecols=(1, 2, 3)
        ).reshape(-1, 3)

    def _read_freqs(idx):
        n = 0
        freqs = []
        stop = 0
        while not freqs or stop < 2:
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

    def _read_electron_phonon(idx):
        results = {}

        broad_re = r"^\s*Gaussian\s*Broadening:\s+([\d.]+)\s+Ry, ngauss=\s+\d+"
        dos_re = r"^\s*DOS\s*=\s*([\d.]+)\s*states/spin/Ry/Unit\s*Cell\s*at\s*Ef=\s+([\d.]+)\s+eV"
        lg_re = r"^\s*lambda\(\s+(\d+)\)=\s+([\d.]+)\s+gamma=\s+([\d.]+)\s+GHz"
        end_re = r"^\s*Number\s*of\s*q\s*in\s*the\s*star\s*=\s+(\d+)$"

        lambdas = []
        gammas = []

        current = None

        n = 1
        while idx + n < n_lines:
            line = fdo_lines[idx + n]

            broad_match = re.match(broad_re, line)
            dos_match = re.match(dos_re, line)
            lg_match = re.match(lg_re, line)
            end_match = re.match(end_re, line)

            if broad_match:
                if lambdas:
                    results[current]["lambdas"] = lambdas
                    results[current]["gammas"] = gammas
                    lambdas = []
                    gammas = []
                current = broad_match[1]
                results[current] = {}
            elif dos_match:
                results[current]["dos"] = float(dos_match[1])
                results[current]["fermi"] = float(dos_match[2])
            elif lg_match:
                lambdas.append(float(lg_match[2]))
                gammas.append(float(lg_match[3]))

            if end_match:
                results[current]["lambdas"] = lambdas
                results[current]["gammas"] = gammas
                break

            n += 1

        return results

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
        ELECTRON_PHONON: _read_electron_phonon,
    }

    iblocks = np.append(output[QPOINTS], n_lines)

    for past, future in zip(iblocks[:-1], iblocks[1:]):
        qpoint = _read_qpoints(past)
        results[qpoint] = {}
        for prop in properties:
            p = (past < output[prop]) & (output[prop] < future)
            selected = output[prop][p]
            if len(selected) == 0:
                continue
            if unique[prop]:
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

    return results


def read_espresso_pw(filename: str | Path) -> dict[str, Any]:
    """
    Wrapper around ase.io.espresso.read_espresso_out to read the output of a pw.x
    and directly return the dictionary of properties to quacc.

    Parameters
    ----------
    filename
        The name of the output file.

    Returns
    -------
    dict
        The results dictionary.
    """
    atoms = _read(filename, format="espresso-out")
    return dict(atoms.calc.properties())


def write_espresso_pw(
    filename: str | Path, atoms: Atoms, properties: dict[str, Any] | None, **kwargs
) -> None:
    """
    Wrapper around ase.io.espresso.write_espresso_in to write the input of a pw.x

    Parameters
    ----------
    filename
        The name of the output file.
    atoms
        The atoms object.
    properties
        Additional ASE properties
    **kwargs
        Additional keyword arguments to pass to ase.io.espresso.write_espresso_in

    Returns
    -------
    None
    """
    _write(filename, atoms, format="espresso-in", properties=properties, **kwargs)


read_espresso_dict = {"ph": read_espresso_ph, "pw": read_espresso_pw}

write_espresso_dict = {"ph": write_espresso_ph, "pw": write_espresso_pw}
