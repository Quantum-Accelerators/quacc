from __future__ import annotations

import re
from gzip import GzipFile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms


def get_pseudopotential_info(
    pp_dict: dict[str, Any], atoms: Atoms
) -> tuple[float, float, dict[str, str]]:
    """
    Function that parses the pseudopotentials and cutoffs from a preset file.
    The cutoffs are taken from the largest value of the cutoffs among the elements
    present.

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


def parse_ph_patterns(root_dir: str | Path, prefix: str) -> dict[int, int]:
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

    fds = Path(root_dir).glob(f"_ph0/{prefix}.phsave/patterns.*.xml*")

    patterns = {}

    qpt_num = r"<QPOINT_NUMBER>(\d+)</QPOINT_NUMBER>"
    repr_count = r"<NUMBER_IRR_REP>(\d+)</NUMBER_IRR_REP"

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
