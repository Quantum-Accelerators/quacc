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
        setups = config["inputs"].get("setups")
        if isinstance(setups, dict):
            config["inputs"]["setups"] = {
                k: get_setup_suffix(k, v) for k, v in setups.items()
            }

    return config


def get_setup_suffix(symbol: Any, setup: Any) -> Any:
    """
    Convert a full POTCAR name (e.g. "Cu_pv") to the suffix (e.g. "_pv") that
    ASE expects for an element-keyed setup. Suffixes (e.g. "_pv") and setups
    keyed by atom index (which ASE expects as full POTCAR names) are returned
    unchanged.

    Parameters
    ----------
    symbol
        The key of the setups dictionary (an element symbol or atom index).
    setup
        The setup for that key.

    Returns
    -------
    Any
        The setup suffix.
    """
    if isinstance(symbol, str) and isinstance(setup, str):
        return setup.removeprefix(symbol)
    return setup
