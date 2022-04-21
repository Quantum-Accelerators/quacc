"""
Basic utility functions
"""
from __future__ import annotations

import os
from shutil import copy
from typing import Any, Dict, List

import yaml
from monty.os.path import zpath
from monty.shutil import decompress_file


def copy_decompress(src_files: List[str], dst) -> None:
    """
    Copy and decompress files from src to dst.

    Parameters
    ----------
    src_files
        List of files to copy and decompress.
    dst
        Destination directory.
    """
    for f in src_files:
        z_path = zpath(f)
        z_file = os.path.basename(z_path)
        if os.path.exists(z_path):
            copy(z_path, os.path.join(dst, z_file))
            decompress_file(os.path.join(dst, z_file))


def merge_dicts(
    d1: Dict[str, Any],
    d2: Dict[str, Any],
    remove_none: bool = False,
    remove_false: bool = False,
    auto_lowercase: bool = True,
) -> Dict[str, Any]:
    """
    Merges two dictionaries into a single dictionary. If both dictionaries
    have the same key, the value from the second dictionary will be used. This
    is done in a case-insensitive manner.

    Parameters
    ----------
    d1
        First dictionary.
    d2
        Second dictionary, which has priority.
    remove_none
        If True, all keys with a value of None in the merged dictionary will be removed.
    remove_false
        If True, all keys with a value of False in the merged dictionary will be removed.
    auto_lowercase
        If True, all keys will be turned into lowercase.
    """
    if auto_lowercase:
        d1 = {k.lower(): v for k, v in d1.items()}
        d2 = {k.lower(): v for k, v in d2.items()}
    d_merged = {**d1, **d2}
    if remove_none:
        d_merged = {k: v for k, v in d_merged.items() if v is not None}
    if remove_false:
        d_merged = {k: v for k, v in d_merged.items() if v is not False}
    return d_merged


def remove_dict_empties(d: Dict[str, Any]) -> Dict[str, Any]:
    """
    For a given dictionary, recursively remove all items that are None
    or are empty lists/dicts.

    Parameters
    ----------
    d
        Dictionary to jsonify

    Returns
    -------
    Dict
        jsonify'd dictionary
    """

    if isinstance(d, dict):
        return {
            k: remove_dict_empties(v)
            for k, v in d.items()
            if v is not None and not (isinstance(v, (dict, list)) and len(v) == 0)
        }
    if isinstance(d, list):
        return [remove_dict_empties(v) for v in d]
    return d


def load_yaml_calc(yaml_path: str) -> Dict[str, Any]:
    """
    Loads a YAML file containing ASE VASP calcultor settings.

    Parameters
    ----------
    yaml_path
        Path to the YAML file.

    Returns
    -------
    Dict
        The calculator configuration (i.e. settings).
    """

    _, ext = os.path.splitext(yaml_path)
    if not ext:
        yaml_path += ".yaml"

    if not os.path.exists(yaml_path):
        raise ValueError(f"Cannot find {yaml_path}.")

    # Load YAML file
    with open(yaml_path, "r") as stream:
        config = yaml.safe_load(stream)

    # Inherit arguments from any parent YAML files
    # but do not overwrite those in the child file.
    for config_arg in config.copy():
        if "parent" in config_arg:
            parent_config = load_yaml_calc(
                os.path.join(os.path.dirname(yaml_path), config[config_arg])
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
