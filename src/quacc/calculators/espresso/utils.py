from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

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
        The parameters dictionary
    binary
        The binary to check (default pw)

    Returns
    -------
    None
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
        if not (ldisp or qplot) and lqdir and is_grid and qpts != (0, 0, 0):
            LOGGER.warning(
                "lqdir is set to True but ldisp and qplot are set to False. the band structure will still be computed at each step. Please set lqdir to False"
            )  # Tom to Andrew: should we be intrusive and change the parameter? Or just warn?
            # This warning will occur if someone use the grid_phonon_flow for a unique q-point
            # and forgets to manually set lqdir to False. (which is set to True by default in the flow)
            # This should be mentionned in the documentation of the flow
