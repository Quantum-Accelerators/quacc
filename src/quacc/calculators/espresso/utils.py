from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.utils.dicts import merge_dicts, remove_dict_nones

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

    if "pseudopotentials" in config:
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
        pp_input_data = {"system": {"ecutwfc": wfc_cutoff, "ecutrho": rho_cutoff}}

    input_data = config.get("input_data")
    input_data = merge_dicts(pp_input_data, input_data)

    atoms_info = config.get("atoms", {})

    if "pbc" in atoms_info:
        atoms.set_pbc(atoms_info["pbc"])

    kpts = config.get("kpts")
    kspacing = config.get("kspacing")

    return_dict = {
        "input_data": input_data,
        "pseudopotentials": pseudopotentials,
        "kspacing": kspacing,
        "atoms": atoms,
        "kpts": kpts,
    }
    return_dict = remove_dict_nones(return_dict)

    if return_dict.get("kpts") == "gamma":
        return_dict[kpts] = None

    return return_dict
