from __future__ import annotations

import warnings

from ase.io.espresso import Namelist


# Stolen from ase.io.espresso and modified to add a parameter 'keys',
#  Maybe there is a better way of doing this?
def construct_namelist(parameters=None, keys=None, warn=False, **kwargs):
    """
    Construct an ordered Namelist containing all the parameters given (as
    a dictionary or kwargs). Keys will be inserted into their appropriate
    section in the namelist and the dictionary may contain flat and nested
    structures. Any kwargs that match input keys will be incorporated into
    their correct section. All matches are case-insensitive, and returned
    Namelist object is a case-insensitive dict.

    If a key is not known to ase, but in a section within `parameters`,
    it will be assumed that it was put there on purpose and included
    in the output namelist. Anything not in a section will be ignored (set
    `warn` to True to see ignored keys).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    The priority of the keys is:
        kwargs[key] > parameters[key] > parameters[section][key]
    Only the highest priority item will be included.

    Parameters
    ----------
    parameters: dict
        Flat or nested set of input parameters.
    warn: bool
        Enable warnings for unused keys.

    Returns
    -------
    input_namelist: Namelist
        pw.x compatible namelist of input parameters.

    """

    # Convert everything to Namelist early to make case-insensitive
    if parameters is None:
        parameters = Namelist()
    else:
        # Maximum one level of nested dict
        # Don't modify in place
        parameters_namelist = Namelist()
        for key, value in parameters.items():
            if isinstance(value, dict):
                parameters_namelist[key] = Namelist(value)
            else:
                parameters_namelist[key] = value
        parameters = parameters_namelist

    # Just a dict
    kwargs = Namelist(kwargs)

    # Final parameter set
    input_namelist = Namelist()
    # Collect
    for section in keys:
        sec_list = Namelist()
        for key in keys[section]:
            # Check all three separately and pop them all so that
            # we can check for missing values later
            if key in parameters.get(section, {}):
                sec_list[key] = parameters[section].pop(key)
            if key in parameters:
                sec_list[key] = parameters.pop(key)
            if key in kwargs:
                sec_list[key] = kwargs.pop(key)
            # Check if there is a key(i) version (no extra parsing)
            for arg_key in list(parameters.get(section, {})):
                if arg_key.split("(")[0].strip().lower() == key.lower():
                    sec_list[arg_key] = parameters[section].pop(arg_key)
            cp_parameters = parameters.copy()
            for arg_key in cp_parameters:
                if arg_key.split("(")[0].strip().lower() == key.lower():
                    sec_list[arg_key] = parameters.pop(arg_key)
            cp_kwargs = kwargs.copy()
            for arg_key in cp_kwargs:
                if arg_key.split("(")[0].strip().lower() == key.lower():
                    sec_list[arg_key] = kwargs.pop(arg_key)

        # Add to output
        input_namelist[section] = sec_list

    unused_keys = list(kwargs)
    # pass anything else already in a section
    for key, value in parameters.items():
        if key in keys and isinstance(value, dict):
            input_namelist[key].update(value)
        elif isinstance(value, dict):
            unused_keys.extend(list(value))
        else:
            unused_keys.append(key)

    if warn and unused_keys:
        warnings.warn("Unused keys: {}".format(", ".join(unused_keys)))

    return input_namelist


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


def namelist_to_string(parameters):
    # This function aim to write basic input keywords to a file
    # For additional cards (binary specifics) this is down in the child class
    # Stolen from ase.io.espresso,
    # maybe this could be made to a function "write_namelist"
    # and imported here? See with the ASE gods.
    pwi = []
    for section in parameters:
        pwi.append(f"&{section.upper()}\n")
        for key, value in parameters[section].items():
            if value is True:
                pwi.append(f"   {key:16} = .true.\n")
            elif value is False:
                pwi.append(f"   {key:16} = .false.\n")
            else:
                pwi.append(f"   {key:16} = {value!r}\n")
        pwi.append("/\n")  # terminate section
    pwi.append("\n")
    return pwi
