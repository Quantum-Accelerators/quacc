"""
Utility functions for file and path handling
"""
from __future__ import annotations

import os
from datetime import datetime
from random import randint
from shutil import copy

import yaml
from monty.io import zopen
from monty.os.path import zpath
from monty.shutil import decompress_file


def check_logfile(logfile: str, check_str: str) -> bool:
    """
    Check if a logfile has a given string (case-insensitive).

    Parameters
    ----------
    logfile : str
        Path to the logfile.
    check_str : str
        String to check for.

    Returns
    -------
    bool
        True if the string is found in the logfile, False otherwise.
    """
    zlog = zpath(logfile)
    with zopen(zlog, "r") as f:
        for line in f:
            if not isinstance(line, str):
                line = line.decode("utf-8")
            if check_str.lower() in line.lower():
                return True
    return False


def copy_decompress(source_files: list[str], destination: str) -> None:
    """
    Copy and decompress files from source to destination.

    Parameters
    ----------
    source_files
        List of files to copy and decompress.
    destination
        Destination directory.

    Returns
    -------
    None
    """
    for f in source_files:
        z_path = zpath(f)
        if os.path.exists(z_path):
            z_file = os.path.basename(z_path)
            copy(z_path, os.path.join(destination, z_file))
            decompress_file(os.path.join(destination, z_file))


def make_job_dir(base_path: str = None) -> str:
    """
    Make a job directory with a unique name. Uses the same format as Jobflow.

    Parameters
    ----------
    base_path
        Path to the base directory. If None, the current working directory is used.

    Returns
    -------
    str
        Path to the job directory.
    """
    if base_path is None:
        base_path = os.getcwd()
    time_now = datetime.utcnow().strftime("%Y-%m-%d-%H-%M-%S-%f")
    job_dir = os.path.join(base_path, f"quacc_{time_now}-{randint(10000, 99999)}")
    os.mkdir(job_dir)

    return job_dir


def load_yaml_calc(yaml_path: str) -> dict:
    """
    Loads a YAML file containing ASE VASP calcultor settings.

    Parameters
    ----------
    yaml_path
        Path to the YAML file.

    Returns
    -------
    dict
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
