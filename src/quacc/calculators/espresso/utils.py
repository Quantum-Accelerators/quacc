from __future__ import annotations

import re
from gzip import GzipFile
from pathlib import Path
from typing import TYPE_CHECKING

from quacc.atoms.core import copy_atoms
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms


def parse_pw_preset(config: dict[str, Any], atoms: Atoms) -> dict[str, Any]:
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
    dict
        A dictionary containing the pseudopotentials and cutoffs
    """

    atoms_copy = copy_atoms(atoms)

    if "pseudopotentials" in config:
        pp_dict = config["pseudopotentials"]
        unique_elements = list(set(atoms_copy.symbols))
        wfc_cutoff, rho_cutoff = 0, 0
        pseudopotentials = {}
        for element in unique_elements:
            if pp_dict[element]["cutoff_wfc"] > wfc_cutoff:
                wfc_cutoff = pp_dict[element]["cutoff_wfc"]
            if pp_dict[element]["cutoff_rho"] > rho_cutoff:
                rho_cutoff = pp_dict[element]["cutoff_rho"]
            pseudopotentials[element] = pp_dict[element]["filename"]
        pp_input_data = {"system": {"ecutwfc": wfc_cutoff, "ecutrho": rho_cutoff}}

    input_data = config.get("input_data")
    input_data = recursive_dict_merge(pp_input_data, input_data)

    kspacing = config.get("kspacing")

    kpts = config.get("kpts")

    return {
        "input_data": input_data,
        "pseudopotentials": pseudopotentials,
        "kspacing": kspacing,
        "kpts": kpts,
    }

  
def parse_ph_patterns(root_dir):
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
