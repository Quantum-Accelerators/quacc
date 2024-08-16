"""I/O utilities for the Vasp calculator."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any


def load_vasp_yaml_calc(yaml_path: str | Path) -> dict[str, Any]:
    """
    Loads a YAML file containing calculator settings. Used for VASP calculations
    and can read quacc-formatted YAMLs that are of the following format:
    ```yaml
    inputs:
      xc: pbe
      algo: all
      setups:
        Cu: Cu_pv
      elemental_magmoms:
        Fe: 5
        Cu: 1
    ```
    where `inputs` is a dictionary of ASE-style input parameters, `setups` is a
    dictionary of ASE-style pseudopotentials, and and `elemental_magmoms` is a
    dictionary of element-wise initial magmoms.

    Parameters
    ----------
    yaml_path
        Path to the YAML file.

    Returns
    -------
    dict
        The calculator configuration (i.e. settings).
    """
    config = load_yaml_calc(yaml_path)

    # Allow for either "Cu_pv" and "_pv" style setups
    if "inputs" in config:
        config["inputs"] = {
            k.lower(): v.lower() if isinstance(v, str) else v
            for k, v in config["inputs"].items()
        }
        for k, v in config["inputs"].get("setups", {}).items():
            if k in v:
                config["inputs"]["setups"][k] = v.split(k)[-1]

    return config
