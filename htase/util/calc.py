import yaml
import os
import numpy as np
import os
from copy import deepcopy


def load_yaml_calc(yaml_file, default_calcs_dir=None):
    """
    Loads a YAML file containing ASE VASP calcultor settings.

    Args:
        yaml_file (str): Filename or full path to YAML file.
        default_calcs_dir (str): If yaml_file is just a filename, load_yaml_calc
            will look in default_calcs_dir for the file.

    Returns:
        calc_preset (dict): The calculator configuration (i.e. settings).
    """

    _, ext = os.path.splitext(yaml_file)
    if not ext:
        yaml_file += ".yaml"

    if os.path.exists(yaml_file):
        yaml_path = yaml_file
    elif default_calcs_dir and os.path.exists(
        os.path.join(default_calcs_dir, yaml_file)
    ):
        yaml_path = os.path.join(default_calcs_dir, yaml_file)
    else:
        raise ValueError(
            f"Cannot find {yaml_file}. Provide the full path or place it in {default_calcs_dir}."
        )

    # Load YAML file
    with open(yaml_path, "r") as stream:
        config = yaml.safe_load(stream)

    # Inherit arguments from any parent YAML files
    # but do not overwrite those in the child file.
    parent_args = ["parent", "parent_magmoms", "parent_setups"]
    for config_arg in parent_args:
        if config_arg in config:
            parent_config = load_yaml_calc(
                config[config_arg], default_calcs_dir=default_calcs_dir
            )
            for k, v in parent_config.items():
                if k not in config:
                    config[k] = v
                else:
                    v_new = parent_config.get(k, {})
                    for kk, vv in v_new.items():
                        if kk not in config[k]:
                            config[k][kk] = vv

    # Allow for either "Cu_pv" and "_pv" style setups
    for k, v in config["inputs"].get("setups", {}).items():
        if k in v:
            config["inputs"]["setups"][k] = v.split(k)[-1]

    return config


def set_magmoms(
    atoms, elemental_mags_dict=None, copy_magmoms=True, mag_default=1.0, mag_cutoff=0.05
):
    """
    Sets the initial magnetic moments in the Atoms object.

    This function deserves particular attention. The following logic is applied:
    - If there is a converged set of magnetic moments, those are moved to the
    initial magmoms if copy_magmoms is True.
    - If there is no converged set of magnetic moments but the user has set initial magmoms,
    those are simply used as is.
    - If there are no converged magnetic moments or initial magnetic moments, then
    the default magnetic moments from the preset (if specified) are set as the
    initial magnetic moments.
    - For any of the above scenarios, if mag_cutoff is not None, the newly set
    initial magnetic moments are checked. If all have a magnitude below mag_cutoff,
    then they are all set to 0 (no spin polarization).

    Args:
        atoms (ase.Atoms): Atoms object
        elemental_mags_dict (dict): Dictionary of elements and their
            corresponding magnetic moments to set.
            Default: None.
        copy_magmoms (bool): Whether to copy the magnetic moments from the
            converged set of magnetic moments to the initial magnetic moments.
            Default: True.
        mag_default (float): Default magnetic moment to use if no magnetic
            moments are specified in the preset.
            Default: 1.0.
        mag_cutoff (float): Magnitude below which the magnetic moments are
            considered to be zero.
            Default: 0.05.

    Returns:
        atoms (ase.Atoms): Atoms object

    """
    atoms = deepcopy(atoms)

    # Handle the magnetic moments
    # Check if a prior job was run and pull the prior magmoms
    if hasattr(atoms, "calc") and getattr(atoms.calc, "results", None) is not None:
        mags = atoms.calc.results.get("magmoms", [0.0] * len(atoms))
        # Note: It is important that we set mags to 0.0 here rather than None if the
        # calculator has no magmoms because: 1) ispin=1 might be set, and 2) we do
        # not want the preset magmoms to be used.
    else:
        mags = None

    # Check if the user has set any initial magmoms
    has_initial_mags = atoms.has("initial_magmoms")

    # If there are no initial magmoms set and this is not a follow-up job,
    # we may need to add some from the preset yaml.
    if mags is None:
        if not has_initial_mags:

            # If the preset dictionary has default magmoms, set
            # those by element. If the element isn't in the magmoms dict
            # then set it to mag_default.
            if elemental_mags_dict:
                initial_mags = np.array(
                    [
                        elemental_mags_dict.get(atom.symbol, mag_default)
                        for atom in atoms
                    ]
                )
                atoms.set_initial_magnetic_moments(initial_mags)
    # Copy converged magmoms to input magmoms, if copy_magmoms is True
    else:
        if copy_magmoms:
            atoms.set_initial_magnetic_moments(mags)

    # If all the set mags are below mag_cutoff, set them to 0
    if mag_cutoff:
        has_new_initial_mags = atoms.has("initial_magmoms")
        new_initial_mags = atoms.get_initial_magnetic_moments()
        if has_new_initial_mags and np.all(np.abs(new_initial_mags) < mag_cutoff):
            atoms.set_initial_magnetic_moments([0.0] * len(atoms))

    return atoms
