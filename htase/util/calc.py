import yaml
import os


def cache_calc(atoms):
    """
    Store the calculator results in atoms.info["results"] for later retrieval.
    This makes it so that things like the converged magnetic moments are not
    lost between serialize/deserialize cycles and also makes it possible to
    retain information about computed properties when manipulating the
    Atoms object, e.g. after a supercell transformation.

    cache_calc supports multiple prior calculators. Each one will be stored in
    atoms.info["results"] = {"calc0": {}, "calc1": {}, ...} with higher numbers
    being the most recent.

    Args:
        atoms (ase.Atoms): Atoms object

    Returns:
        atoms (ase.Atoms): Atoms object with calculator results attached in atoms.info["results"]
    """
    if hasattr(atoms, "calc") and hasattr(atoms.calc, "results"):
        if atoms.info.get("results", None) is None:
            prior_calcs = 0
        else:
            prior_calcs = len(atoms.info["results"])

        if atoms.info.get("results", None) is None:
            atoms.info["results"] = {}
        atoms.info["results"][f"calc{prior_calcs}"] = atoms.calc.results

    return atoms


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
    for k, v in config["inputs"].get("setups", {}).items():
        if k in v:
            config["inputs"]["setups"][k] = v.split(k)[-1]

    return config
