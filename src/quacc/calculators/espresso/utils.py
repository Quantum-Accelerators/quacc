from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.io.espresso import Namelist

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
