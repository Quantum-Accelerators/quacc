from __future__ import annotations

import re
from logging import getLogger
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from quacc.utils.dicts import Remove

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from ase.atoms import Atoms
    from numpy.typing import ArrayLike

    from quacc.types import SourceDirectory

LOGGER = getLogger(__name__)


def matdyn_qpoint_cards(qpts: ArrayLike, q_in_band_form: bool = False) -> list[str]:
    """Format matdyn q-points for ASE's ``additional_cards`` writer.

    Explicit points have three columns. Band-form points have a fourth column
    containing QE's interpolation count: ``n`` intervals to the next endpoint,
    or zero for a discontinuity. QE ignores the last endpoint's count.
    """
    points = np.asarray(qpts, dtype=float)
    ncols = 4 if q_in_band_form else 3
    if points.ndim != 2 or points.shape[1] != ncols or len(points) == 0:
        raise ValueError(f"matdyn qpts must be a nonempty (n, {ncols}) array.")
    if not np.isfinite(points).all():
        raise ValueError("matdyn qpts must contain only finite values.")
    if q_in_band_form:
        counts = points[:, 3]
        if np.any(counts < 0) or np.any(counts != np.floor(counts)):
            raise ValueError(
                "matdyn interpolation counts must be nonnegative integers."
            )

    cards = [str(len(points))]
    for point in points:
        line = " ".join(f"{value:.16g}" for value in point[:3])
        if q_in_band_form:
            line += f" {int(point[3])}"
        cards.append(line)
    return cards


def read_ifc_has_zstar(filename: str | Path) -> bool | None:
    """Read whether q2r.x force constants contain the dielectric tensor and Z*.

    This is QE's ``has_zstar`` flag. q2r.x sets it when the Gamma dynamical
    matrix includes the dielectric tensor and Born effective charges, which
    ph.x computes for insulators. matdyn.x then adds the dipole-dipole term
    and, at Gamma, the non-analytic (LO-TO) term. Only text files are read.
    None is returned for XML, missing, or unrecognized files: matdyn.x has
    already read the file, so its results remain valid without this flag.
    """
    path = Path(filename)
    if path.suffix == ".xml" or not path.is_file():
        return None
    with path.open() as fd:
        try:
            # The first record starts with ntyp, nat, ibrav, and celldm(1:6).
            ntyp, nat, ibrav = (int(value) for value in fd.readline().split()[:3])
            # Skip the lattice vectors (ibrav=0 only), species, and positions.
            for _ in range(ntyp + nat + (3 if ibrav == 0 else 0)):
                fd.readline()
            flag = fd.readline().split()[0].strip(".").upper()[:1]
        except (ValueError, IndexError):
            flag = ""
    if flag not in ("T", "F"):
        LOGGER.warning(f"Could not read has_zstar from the header of {filename}.")
        return None
    return flag == "T"


def read_matdyn_freq(filename: str | Path) -> dict[str, Any]:
    """Read matdyn frequencies and, if present, the companion ``.gp`` table.

    Frequencies retain QE's mode and q-point ordering and signs, in cm^-1.
    The ``&plot`` file stores Cartesian q-points in units of 2pi/alat, even
    when the input is fractional. Distances are read only from ``.gp``; QE's
    cumulative distances include jumps between disconnected path segments.
    If only ``.gp`` exists, q-point coordinates are unavailable and omitted.
    """
    filename = Path(filename)
    results: dict[str, Any] = {}
    if filename.is_file():
        content = filename.read_text()
        # ASE's namelist reader includes the section token in the first inline
        # assignment key ("&plot nbnd"), so read QE's compact header directly.
        header = re.match(r"\s*&plot\b([^/]+)/", content, re.IGNORECASE)
        if header is None:
            raise ValueError(f"Invalid matdyn frequency header in {filename}.")
        dimensions = {
            key.lower(): int(value)
            for key, value in re.findall(
                r"\b(nbnd|nks)\s*=\s*([+-]?\d+)\b", header[1], re.IGNORECASE
            )
        }
        try:
            nbnd, nks = dimensions["nbnd"], dimensions["nks"]
        except KeyError as exc:
            raise ValueError(f"Invalid matdyn frequency header in {filename}.") from exc
        if nbnd <= 0 or nks <= 0:
            raise ValueError(f"Invalid matdyn frequency dimensions in {filename}.")

        qpoints = np.empty((nks, 3))
        frequencies = np.empty((nks, nbnd))
        records = iter(
            line for line in content[header.end() :].splitlines() if line.strip()
        )
        try:
            for iq in range(nks):
                point = _matdyn_numbers(next(records))
                # QE versions write either three coordinates or an additional
                # weight/path coordinate. The latter is not a fractional q-point.
                if len(point) not in (3, 4):
                    raise ValueError("Expected three matdyn q-point coordinates.")
                qpoints[iq] = point[:3]
                modes: list[float] = []
                while len(modes) < nbnd:
                    modes.extend(_matdyn_numbers(next(records)))
                if len(modes) != nbnd:
                    raise ValueError("Incorrect number of matdyn frequencies.")
                frequencies[iq] = modes
        except StopIteration as exc:
            raise ValueError(
                f"Incomplete matdyn frequency data in {filename}."
            ) from exc
        if any(line.strip() for line in records):
            raise ValueError(f"Unexpected extra matdyn frequency data in {filename}.")
        results.update(
            qpoints=qpoints,
            qpoint_coordinate_system="cartesian",
            qpoint_units="2pi/alat",
            frequencies=frequencies,
        )

    gp_file = filename.with_name(filename.name + ".gp")
    if gp_file.is_file():
        rows = [
            _matdyn_numbers(line)
            for line in gp_file.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        try:
            table = np.asarray(rows, dtype=float)
        except ValueError as exc:
            raise ValueError(f"Inconsistent matdyn .gp rows in {gp_file}.") from exc
        if table.ndim != 2 or table.shape[1] < 2:
            raise ValueError(f"Invalid matdyn .gp table in {gp_file}.")
        if "frequencies" in results and (
            table[:, 1:].shape != results["frequencies"].shape
            or not np.allclose(table[:, 1:], results["frequencies"], rtol=0, atol=1e-4)
        ):
            raise ValueError(f"matdyn frequency and .gp files disagree: {filename}.")
        results.setdefault("frequencies", table[:, 1:])
        results.update(distances=table[:, 0], distance_units="2pi/alat")

    if results:
        results["frequency_units"] = "cm^-1"
    return results


def _matdyn_numbers(line: str) -> list[float]:
    """Read QE numeric records, including Fortran exponents and adjacent signs."""
    pattern = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
    if re.sub(pattern, "", line).strip():
        raise ValueError(f"Invalid numeric matdyn record: {line!r}")
    return [
        float(value.replace("D", "e").replace("d", "e"))
        for value in re.findall(pattern, line)
    ]


def get_pseudopotential_info(
    pp_dict: dict[str, Any], atoms: Atoms
) -> tuple[float, float, dict[str, str]]:
    """
    Function that parses the pseudopotentials and cutoffs from a preset file. The
    cutoffs are taken from the largest value of the cutoffs among the elements present.

    Parameters
    ----------
    pp_dict
        The "pseudopotential" parameters returned by load_yaml_calc
    atoms
        The atoms object to get the elements from

    Returns
    -------
    float
        The max(ecutwfc) value
    float
        The max(ecutrho) value
    dict[str, str]
        The pseudopotentials dictionary, e.g. {"O": "O.pbe-n-kjpaw_psl.0.1.UPF"}
    """
    unique_elements = list(set(atoms.get_chemical_symbols()))
    ecutwfc, ecutrho = 0, 0
    pseudopotentials = {}
    for element in unique_elements:
        if pp_dict[element]["cutoff_wfc"] > ecutwfc:
            ecutwfc = pp_dict[element]["cutoff_wfc"]
        if pp_dict[element]["cutoff_rho"] > ecutrho:
            ecutrho = pp_dict[element]["cutoff_rho"]
        pseudopotentials[element] = pp_dict[element]["filename"]
    return ecutwfc, ecutrho, pseudopotentials


def grid_copy_files(
    ph_input_data: Mapping[str, Any],
    directory: str | Path,
    qnum: int | str,
    qpt: tuple[float, float, float],
) -> dict[SourceDirectory, list[Path]]:
    """
    Function that returns a dictionary of files to copy for the grid calculation.

    Parameters
    ----------
    ph_input_data
        The input data for the ph calculation
    directory
        The directory name to copy the files from
    qnum
        The q-point number
    qpt
        The q-point coordinates in QE units

    Returns
    -------
    dict[SourceDirectory, Filenames]
        The dictionary of files to copy
    """
    lqdir = ph_input_data["inputph"].get("lqdir", False)

    files_to_copy = {
        directory: [
            Path("_ph0", "pwscf.phsave", "control_ph.xml*"),
            Path("_ph0", "pwscf.phsave", "status_run.xml*"),
            Path("_ph0", "pwscf.phsave", "patterns.*.xml*"),
            Path("_ph0", "pwscf.phsave", "tensors.xml*"),
        ]
    }

    if lqdir or qpt == (0.0, 0.0, 0.0):
        files_to_copy[directory].extend(
            [
                Path("pwscf.save", "charge-density.*"),
                Path("pwscf.save", "data-file-schema.xml*"),
                Path("pwscf.save", "paw.txt*"),
                Path("pwscf.save", "wfc*.*"),
            ]
        )
        if qpt != (0.0, 0.0, 0.0):
            files_to_copy[directory].extend(
                [
                    Path("_ph0", f"pwscf.q_{qnum}", "pwscf.save", "*"),
                    Path("_ph0", f"pwscf.q_{qnum}", "pwscf.wfc*"),
                ]
            )
    else:
        files_to_copy[directory].extend(
            [Path("_ph0", "pwscf.wfc*"), Path("_ph0", "pwscf.save", "*")]
        )

    return files_to_copy


def grid_prepare_repr(patterns: dict[str, Any], nblocks: int) -> list:
    """
    Function that prepares the representations for the grid calculation.

    Parameters
    ----------
    patterns
        The representation patterns dictionary from the ph_init_job_results
    nblocks
        The number of blocks to group the representations in. See
        [quacc.recipes.espresso.phonons.grid_phonon_flow][].

    Returns
    -------
    list
        The list of representations to do grouped in blocks if nblocks > 1
    """
    this_block = nblocks if nblocks > 0 else len(patterns)
    repr_to_do = [rep for rep in patterns if not patterns[rep]["done"]]
    return np.array_split(repr_to_do, np.ceil(len(repr_to_do) / this_block))


def espresso_prepare_dir(outdir: str | Path, binary: str = "pw") -> dict[str, Any]:
    """
    Function that prepares the espresso dictionary for the calculation.

    Parameters
    ----------
    outdir
        The output to be used for the espresso calculation
    binary
        The binary to be used for the espresso calculation

    Returns
    -------
    dict[str, Any]
        Input data for the espresso calculation
    """
    outkeys: dict[str, dict[str, Any]] = {
        "pw": {"control": {"prefix": "pwscf", "outdir": outdir, "wfcdir": Remove}},
        "ph": {
            "inputph": {
                "prefix": "pwscf",
                "fildyn": "matdyn",
                "outdir": outdir,
                "ahc_dir": Remove,
                "wpot_dir": Remove,
                "dvscf_star%dir": Remove,
                "drho_star%dir": Remove,
            }
        },
        "pp": {"inputpp": {"prefix": "pwscf", "filplot": "tmp.pp", "outdir": outdir}},
        "dos": {"dos": {"prefix": "pwscf", "fildos": "pwscf.dos", "outdir": outdir}},
        "projwfc": {
            "projwfc": {"prefix": "pwscf", "filpdos": "pwscf", "outdir": outdir}
        },
        "matdyn": {
            "input": {
                "flfrc": "q2r.fc",
                "fldos": "matdyn.dos",
                "flfrq": "matdyn.freq",
                "flvec": "matdyn.modes",
                "fleig": "matdyn.eig",
            }
        },
        "q2r": {"input": {"fildyn": "matdyn", "flfrc": "q2r.fc"}},
        "bands": {
            "bands": {"prefix": "pwscf", "filband": "bands.dat", "outdir": outdir}
        },
        "fs": {
            "fermi": {
                "prefix": "pwscf",
                "file_fs": "fermi_surface.bxsf",
                "outdir": outdir,
            }
        },
        "dvscf_q2r": {
            "input": {
                "prefix": "pwscf",
                "fildyn": "matdyn",
                "outdir": outdir,
                "wpot_dir": Remove,
            }
        },
        "postahc": {"input": {"ahc_dir": "ahc_dir/", "flvec": "matdyn.modes"}},
    }

    return outkeys.get(binary, {})


def prepare_copy_files(parameters: dict[str, Any], binary: str = "pw") -> list[Path]:
    """
    Function that prepares the copy files for the espresso calculation.

    Parameters
    ----------
    parameters
        The input data for the espresso calculation
    binary
        The binary to use for the espresso calculation

    Returns
    -------
    list[Path]
        Paths to copy for the espresso calculation
    """
    to_copy = []

    pw_base = [
        Path("pwscf.save", "charge-density.*"),
        Path("pwscf.save", "data-file-schema.*"),
        Path("pwscf.save", "paw.*"),
    ]

    input_data = parameters.get("input_data", {})

    if binary == "pw":
        control = input_data.get("control", {})
        restart_mode = control.get("restart_mode", "from_scratch")
        electrons = input_data.get("electrons", {})
        startingpot = electrons.get("startingpot", "atomic")
        startingwfc = electrons.get("startingwfc", "atomic+random")
        calculation = control.get("calculation", "scf")

        if restart_mode == "restart":
            to_copy.extend(
                [
                    Path("pwscf.wfc*"),
                    Path("pwscf.mix*"),
                    Path("pwscf.restart_k*"),
                    Path("pwscf.restart_scf*"),
                ]
            )

        need_chg_dens = (
            startingpot == "file"
            or calculation in ["bands", "nscf"]
            or restart_mode == "restart"
        )

        need_wfc = startingwfc == "file" or restart_mode == "restart"

        if need_chg_dens:
            to_copy.append(Path("pwscf.save", "charge-density.*"))

        if need_wfc:
            to_copy.append(Path("pwscf.save", "wfc*.*"))

        to_copy.extend(
            [Path("pwscf.save", "data-file-schema.*"), Path("pwscf.save", "paw.*")]
        )

    elif binary in {"ph", "phcg"}:
        to_copy.extend(pw_base)
        to_copy.append(Path("pwscf.save", "wfc*.*"))

        inputph = input_data.get("inputph", {})
        ldisp = inputph.get("ldisp", False)
        fildvscf = inputph.get("fildvscf", "")
        recover = inputph.get("recover", False)
        lqdir = inputph.get("lqdir", False) or (ldisp and fildvscf)
        ldvscf_interpolate = inputph.get("ldvscf_interpolate", False)

        if lqdir:
            to_copy.extend(
                [Path("_ph*", "pwscf.q_*", "pwscf.save", "data-file-schema.*")]
            )

            if recover:
                to_copy.extend(
                    [Path("_ph*", "pwscf.q_*", "pwscf.save", "charge-density.*")]
                )

        if recover:
            to_copy.append(Path("_ph*", "pwscf.phsave"))

        if ldvscf_interpolate:
            to_copy.extend((Path("_ph*", "pwscf.dvscf*"), Path("w_pot")))
            if lqdir:
                to_copy.append(Path("_ph*", "pwscf.q_*", "pwscf.dvscf*"))

    elif binary in {"dos", "fs"}:
        to_copy.extend(pw_base)

    elif binary in {"projwfc", "bands"}:
        to_copy.extend(pw_base)
        to_copy.append(Path("pwscf.save", "wfc*.*"))

    elif binary == "pp":
        plotnum = input_data.get("plot_num", 0)
        wfc_needed = [3, 7, 10]

        to_copy.extend(pw_base)

        if plotnum in wfc_needed:
            to_copy.append(Path("pwscf.save", "wfc*.*"))

    elif binary == "matdyn":
        to_copy.append(Path("q2r.fc*"))

    elif binary == "q2r":
        to_copy.append(Path("matdyn*"))

    elif binary == "dvscf_q2r":
        to_copy.extend(pw_base)
        to_copy.extend(
            [
                Path("matdyn0*"),
                Path("_ph*", "pwscf.phsave"),
                Path("_ph*", "pwscf.dvscf*"),
                Path("_ph*", "pwscf.q_*", "pwscf.dvscf*"),
            ]
        )

    elif binary == "postahc":
        to_copy.extend(pw_base)
        to_copy.extend([Path("ahc_dir"), Path("matdyn.modes*")])

    return to_copy


def remove_conflicting_kpts_kspacing(
    to_change_dict: dict[str, Any], reference_dict: dict[str, Any]
) -> dict[str, Any]:
    """
    Parameters
    ----------
    to_change_dict
        The dictionary that will be changed
    reference_dict
        The reference dictionary that will remain intact

    Returns
    -------
    dict[str, Any]
        The modified dictionary
    """
    if "kpts" in reference_dict:
        to_change_dict.pop("kspacing", None)
    if "kspacing" in reference_dict:
        to_change_dict.pop("kpts", None)

    return to_change_dict
