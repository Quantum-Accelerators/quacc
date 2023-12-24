from __future__ import annotations

import re
from gzip import GzipFile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms


def parse_pp_and_cutoff(config: dict[str, Any], atoms: Atoms) -> dict[str, Any] | None:
    """
    Function that parses the pseudopotentials and cutoffs from a preset file.
    The cutoffs are taken from the largest value of the cutoffs among the elements
    present.

    Parameters
    ----------
    config
        The config dictionary return by load_yaml_calc
    atoms
        The atoms object to get the elements from

    Returns
    -------
    dict | None
        A dictionary containing the pseudopotentials and cutoffs
    """

    if "pseudopotentials" not in config:
        return None

    pp_dict = config["pseudopotentials"]
    unique_elements = list(set(atoms.symbols))
    wfc_cutoff, rho_cutoff = 0, 0
    pseudopotentials = {}
    for element in unique_elements:
        if pp_dict[element]["cutoff_wfc"] > wfc_cutoff:
            wfc_cutoff = pp_dict[element]["cutoff_wfc"]
        if pp_dict[element]["cutoff_rho"] > rho_cutoff:
            rho_cutoff = pp_dict[element]["cutoff_rho"]
        pseudopotentials[element] = pp_dict[element]["filename"]
    tmp_input_data = {"system": {"ecutwfc": wfc_cutoff, "ecutrho": rho_cutoff}}

    return {"pseudopotentials": pseudopotentials, "input_data": tmp_input_data}


def parse_ph_patterns(root_dir : str | Path) -> dict[int, int]:
    """
    Function that parses the patterns from a ph.x calculation.

    Parameters
    ----------
    root_dir
        The root directory of the calculation

    Returns
    -------
    dict
        A dictionary containing the patterns
    """
    # Patterns (which means number of representation per q-point)
    # are found in files which follow the pattern:
    fds = Path(root_dir).glob("_ph0/*.phsave/patterns.*.xml*")
    # we do not use the xml parser because of security issues
    patterns = {}

    qpt_num = r"<QPOINT_NUMBER>(\d+)</QPOINT_NUMBER>"
    repr_count = r"<NUMBER_IRR_REP>(\d+)</NUMBER_IRR_REP"
    # If compressed, decompress first in memory
    for fd in fds:
        if fd.suffix == ".gz":
            decompressed = GzipFile(filename=fd, mode="r")
            lines = str(decompressed.read())
        else:
            lines = fd.read_text()
        q = int(re.search(qpt_num, lines, re.MULTILINE).group(1))
        c = int(re.search(repr_count, lines, re.MULTILINE).group(1))

        patterns[q] = c

    return patterns
