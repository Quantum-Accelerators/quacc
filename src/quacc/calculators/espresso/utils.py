from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.io.espresso import Namelist

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.utils.files import Filenames, SourceDirectory

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
) -> dict[SourceDirectory, Filenames]:
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

    files_to_copy = {prev_dir: []}

    basics_to_copy = ["charge-density.*", "data-file-schema.*", "paw.*"]

    if restart_mode == "restart":
        files_to_copy[prev_dir].append(Path(wfcdir, f"{prefix}.wfc*"))
        files_to_copy[prev_dir].append(Path(wfcdir, f"{prefix}.mix*"))
        files_to_copy[prev_dir].append(Path(wfcdir, f"{prefix}.restart_k*"))
        files_to_copy[prev_dir].append(Path(wfcdir, f"{prefix}.restart_scf*"))
    elif include_wfc:
        basics_to_copy.append("wfc*.*")

    files_to_copy[prev_dir].extend(
        [Path(outdir, f"{prefix}.save", i) for i in basics_to_copy]
    )

    return files_to_copy


def grid_copy_files(
    ph_input_data: dict[str, Any],
    dir_name: str | Path,
    qnum: int,
    qpt: tuple[float, float, float],
) -> dict[SourceDirectory, Filenames]:
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

    files_to_copy = {
        dir_name: [
            Path(outdir, "_ph0", f"{prefix}.phsave", "control_ph.xml*"),
            Path(outdir, "_ph0", f"{prefix}.phsave", "status_run.xml*"),
            Path(outdir, "_ph0", f"{prefix}.phsave", "patterns.*.xml*"),
            Path(outdir, "_ph0", f"{prefix}.phsave", "tensors.xml*"),
        ]
    }

    if lqdir or qpt == (0.0, 0.0, 0.0):
        files_to_copy[dir_name].extend(
            [
                Path(outdir, f"{prefix}.save", "charge-density.*"),
                Path(outdir, f"{prefix}.save", "data-file-schema.xml.*"),
                Path(outdir, f"{prefix}.save", "paw.txt.*"),
                Path(outdir, f"{prefix}.save", "wfc*.*"),
            ]
        )
        if qpt != (0.0, 0.0, 0.0):
            files_to_copy[dir_name].extend(
                [
                    Path(outdir, "_ph0", f"{prefix}.q_{qnum}", f"{prefix}.save", "*"),
                    Path(outdir, "_ph0", f"{prefix}.q_{qnum}", f"{prefix}.wfc*"),
                ]
            )
    else:
        files_to_copy[dir_name].extend(
            [
                Path(outdir, "_ph0", f"{prefix}.wfc*"),
                Path(outdir, "_ph0", f"{prefix}.save", "*"),
            ]
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
