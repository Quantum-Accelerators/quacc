from __future__ import annotations


def parse_pp_and_cutoff(config, atoms):
    # if atoms are here, we compute the highest ecutwfc and ecutrho

    # pseudopotentials should be the only special case
    # in the sense that it contains additional information
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
    else:
        return {}
    tmp_input_data = {"system": {"ecutwfc": wfc_cutoff, "ecutrho": rho_cutoff}}
    return {"pseudopotentials": pseudopotentials, "input_data": tmp_input_data}
