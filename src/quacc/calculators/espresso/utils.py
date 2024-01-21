from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms


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
                    f"{outdir}/_ph0/{prefix}.q_{qnum}/{prefix}.phsave/*",
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

