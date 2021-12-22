import yaml
import os


def load_yaml_calc(file_path):
    """
    Loads a YAML file containing ASE VASP calcultor settings.

    Args:
        file_path (str): Path to YAML file.

    Returns:
        config (dict): The calculator configuration (i.e. settings).
    """

    # Load YAML file
    with open(file_path, "r") as stream:
        config = yaml.safe_load(stream)

    # Inherit arguments from any parent YAML files
    # but do not overwrite those in the child file.
    parent_args = ["parent", "parent_magmoms", "parent_setups"]
    for config_arg in parent_args:
        if config_arg in config:
            parent_config = load_yaml_calc(
                os.path.join(os.path.dirname(file_path), config[config_arg])
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
    for setup in config["inputs"].get("setups", {}):
        for k, v in setup.items():
            if k in v:
                config["inputs"]["setups"][k] = v.split(k)[-1]

    return config
