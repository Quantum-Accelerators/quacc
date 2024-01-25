from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.io.espresso import Namelist

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

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
        The pseudopotentials dictinoary, e.g. {"O": "O.pbe-n-kjpaw_psl.0.1.UPF"}
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


def pw_copy_files(
    input_data: dict[str, Any], prev_dir: str | Path, include_wfc: bool = True
) -> dict[str, list[str | Path]]:
    """
    Function that take care of copying the correct files from a previous pw.x
    to a current pw.x/bands.x/dos.x... calculation. wfc in collected format
    might be optionally ommited.

    Parameters
    ----------
    input_data
        input_data of the current calculation
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to partially
        copy the tree structure of that directory to the working directory
        of this calculation.
    include_wfc
        Whether to include the wfc files or not, for dos.x and bands.x they are
        not needed for example.

    Returns
    -------
    dict
        Dictionary of files to copy. The key is the directory to copy from and
        the value is a list of files to copy from that directory.
    """
    input_data = Namelist(input_data)
    input_data.to_nested(binary="pw")

    control = input_data.get("control", {})

    prefix = control.get("prefix", "pwscf")
    restart_mode = control.get("restart_mode", "from_scratch")

    outdir = control.get("outdir", ".")
    wfcdir = control.get("wfcdir", outdir)

    file_to_copy = {prev_dir: []}

    basics_to_copy = ["charge-density.*", "data-file-schema.*", "paw.*"]

    if restart_mode == "restart":
        file_to_copy[prev_dir].append(f"{wfcdir}/{prefix}.wfc*")
        file_to_copy[prev_dir].append(f"{wfcdir}/{prefix}.mix*")
        file_to_copy[prev_dir].append(f"{wfcdir}/{prefix}.restart_k*")
        file_to_copy[prev_dir].append(f"{wfcdir}/{prefix}.restart_scf*")
    elif include_wfc:
        basics_to_copy.append("wfc*.*")

    file_to_copy[prev_dir].extend(
        [f"{outdir}/{prefix}.save/{i}" for i in basics_to_copy]
    )

    return file_to_copy


def grid_copy_files(
    ph_input_data: dict[str, Any],
    dir_name: str | Path,
    qnum: int,
    qpt: tuple[float, float, float],
) -> dict[str, list[str]]:
    """
    Function that returns a dictionary of files to copy for the grid calculation.

    Parameters
    ----------
    ph_input_data
        The input data for the ph calculation
    dir_name
        The directory name to copy the files from
    qnum
        The q-point number
    qpt
        The q-point coordinates in QE units

    Returns
    -------
    dict
        The dictionary of files to copy
    """

    prefix = ph_input_data["inputph"].get("prefix", "pwscf")
    outdir = ph_input_data["inputph"].get("outdir", ".")
    lqdir = ph_input_data["inputph"].get("lqdir", False)

    file_to_copy = {
        dir_name: [
            f"{outdir}/_ph0/{prefix}.phsave/control_ph.xml*",
            f"{outdir}/_ph0/{prefix}.phsave/status_run.xml*",
            f"{outdir}/_ph0/{prefix}.phsave/patterns.*.xml*",
            f"{outdir}/_ph0/{prefix}.phsave/tensors.xml*",
        ]
    }

    if lqdir or qpt == (0.0, 0.0, 0.0):
        file_to_copy[dir_name].extend(
            [
                f"{outdir}/{prefix}.save/charge-density.*",
                f"{outdir}/{prefix}.save/data-file-schema.xml.*",
                f"{outdir}/{prefix}.save/paw.txt.*",
                f"{outdir}/{prefix}.save/wfc*.*",
            ]
        )
        if qpt != (0.0, 0.0, 0.0):
            file_to_copy[dir_name].extend(
                [
                    f"{outdir}/_ph0/{prefix}.q_{qnum}/{prefix}.save/*",
                    f"{outdir}/_ph0/{prefix}.q_{qnum}/{prefix}.wfc*",
                ]
            )
    else:
        file_to_copy[dir_name].extend(
            [f"{outdir}/_ph0/{prefix}.wfc*", f"{outdir}/_ph0/{prefix}.save/*"]
        )

    return file_to_copy


def grid_prepare_repr(patterns: dict[str, Any], nblocks: int) -> list:
    """
    Function that prepares the representations for the grid calculation.

    Parameters
    ----------
    patterns
        The representation patterns dictionary from the ph_init_job_results
    nblocks
        The number of blocks to group the representations in (default 1) see
        [quacc.recipes.espresso.phonons.grid_phonon_flow][].

    Returns
    -------
    list
        The list of representations to do grouped in blocks if nblocks > 1
    """

    this_block = nblocks if nblocks > 0 else len(patterns)
    repr_to_do = [rep for rep in patterns if not patterns[rep]["done"]]
    return np.array_split(repr_to_do, np.ceil(len(repr_to_do) / this_block))


def sanity_checks(parameters: dict[str, Any], binary: str = "pw") -> None:
    """Function that performs sanity checks on the input_data. It is meant
    to catch common mistakes that are not caught by the espresso binaries.

    Parameters
    ----------
    parameters
        The parameters dictionary which is assumed to already be in
        the nested format.
    binary
        The binary to check (default pw).

    Returns
    -------
    dict
        The modified parameters dictionary.
    """

    input_data = parameters.get("input_data", {})

    if binary == "ph":
        input_ph = input_data.get("inputph", {})
        qpts = parameters.get("qpts", (0, 0, 0))

        qplot = input_ph.get("qplot", False)
        lqdir = input_ph.get("lqdir", False)
        recover = input_ph.get("recover", False)
        ldisp = input_ph.get("ldisp", False)

        is_grid = input_ph.get("start_q") or input_ph.get("start_irr")
        # Temporary patch for https://gitlab.com/QEF/q-e/-/issues/644
        if qplot and lqdir and recover and is_grid:
            prefix = input_ph.get("prefix", "pwscf")
            outdir = input_ph.get(
                "outdir", "."
            )  # TODO: change to self.outdir when DOS PR is merged
            Path(outdir, "_ph0", f"{prefix}.q_1").mkdir(parents=True, exist_ok=True)
        if not (ldisp or qplot):
            if np.array(qpts).shape == (1, 4):
                LOGGER.warning(
                    "qpts is a 2D array despite ldisp and qplot being set to False. Converting to 1D array"
                )
                qpts = tuple(qpts[0])
            if lqdir and is_grid and qpts != (0, 0, 0):
                LOGGER.warning(
                    "lqdir is set to True but ldisp and qplot are set to False. The band structure will still be computed at each step. Setting lqdir to False"
                )
                input_ph["lqdir"] = False

        parameters["input_data"]["inputph"] = input_ph
        parameters["qpts"] = qpts

    return parameters
