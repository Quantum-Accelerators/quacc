import yaml
import os
from copy import deepcopy


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

    Also moves .get_magnetic_moments() to .get_initial_magnetic_moments()

    Args:
        atoms (ase.Atoms): Atoms object

    Returns:
        atoms (ase.Atoms): Atoms object with calculator results attached in atoms.info["results"]
    """
    atoms = deepcopy(atoms)

    if hasattr(atoms, "calc") and getattr(atoms.calc, "results", None) is not None:

        # Dump calculator results into the .info tag
        atoms.calc.results["rundir"] = os.getcwd()
        if atoms.info.get("results", None) is None:
            prior_calcs = 0
            atoms.info["results"] = {}
        else:
            prior_calcs = len(atoms.info["results"])

        atoms.info["results"][f"calc{prior_calcs}"] = atoms.calc.results

        # Move converged magmoms to initial magmoms
        # If none were present, then initial magmoms should be set to 0's
        atoms.set_initial_magnetic_moments(
            atoms.calc.results.get("magmoms", [0.0] * len(atoms))
        )

        # Clear off the calculator so we can run a new job
        atoms.calc = None

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
