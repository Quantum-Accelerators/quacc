from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from quacc.utils.dicts import Remove

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, SourceDirectory

LOGGER = logging.getLogger(__name__)


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
    ph_input_data: dict[str, Any],
    directory: str | Path,
    qnum: int,
    qpt: tuple[float, float, float],
) -> dict[SourceDirectory, Filenames]:
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
    outkeys = {
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
            "bands": {"prefix": "pwscf", "filband": "bands.out", "outdir": outdir}
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
