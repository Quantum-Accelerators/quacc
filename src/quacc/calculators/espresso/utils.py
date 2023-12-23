from __future__ import annotations

from typing import TYPE_CHECKING

from ase.io.espresso import kspacing_to_grid

from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms


def parse_pw_preset(config: dict[str, Any], atoms: Atoms) -> dict[str, Any] | None:
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

    atoms_copy = atoms.copy()

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
