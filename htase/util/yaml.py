from typing import Optional, Dict
import yaml
import os


def load_yaml_calc(yaml_file: str, default_calcs_dir: Optional[str] = None) -> Dict:
    """
    Loads a YAML file containing ASE VASP calcultor settings.

    Parameters
    ----------
    yaml_file
        Filename or full path to YAML file.
    default_calcs_dir
        If yaml_file is just a filename, load_yaml_calc will look in default_calcs_dir for the file.

    Returns
    -------
    Dict
        The calculator configuration (i.e. settings).
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


def load_yaml_settings(yaml_file: str) -> Dict:
    """
    Loads a standard YAML settings file. Any entry marked with
    a "$" sign is an environment variable.

    Parameters
    ----------
    yaml_file
        Filename or full path to YAML file.

    Returns
    -------
    Dict
        The settings specified in the YAML file.
    """

    settings = yaml.safe_load(open(yaml_file))

    # If $ is the first character, get from the environment variable
    for k, v in settings.items():
        if isinstance(v, str) and v[0] == "$":
            if os.path.expandvars(v) == v:
                raise EnvironmentError(
                    f"Missing environment variable {v}, as specified in {yaml_file}"
                )
            settings[k] = os.path.expandvars(v)

    return settings
