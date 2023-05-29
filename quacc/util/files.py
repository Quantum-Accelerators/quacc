"""
Basic utility functions
"""
from __future__ import annotations

import contextlib
import os
import socket
from pathlib import Path
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


def find_recent_logfile(dir_name: str, logfile_extensions: str | list[str]):
    """
    Find the most recent logfile in a given directory.

    Parameters
    ----------
    dir_name
        The path to the directory to search
    logfile_extensions
        The extension (or list of possible extensions) of the logfile to search for.
        For an exact match only, put in the full file name.

    Returns
    -------
    logfile
        The path to the most recent logfile with the desired extension
    """
    mod_time = 0.0
    logfile = None
    if isinstance(logfile_extensions, str):
        logfile_extensions = [logfile_extensions]
    for f in os.listdir(dir_name):
        f_path = os.path.join(dir_name, f)
        for ext in logfile_extensions:
            if ext in f and os.path.getmtime(f_path) > mod_time:
                mod_time = os.path.getmtime(f_path)
                logfile = os.path.abspath(f_path)
    return logfile


def copy_decompress(src_files: list[str], dst: str) -> None:
    """
    Copy and decompress files from src to dst.

    Parameters
    ----------
    src_files
        List of files to copy and decompress.
    dst
        Destination directory.

    Returns
    -------
    None
    """
    for f in src_files:
        z_path = zpath(f)
        if os.path.exists(z_path):
            z_file = os.path.basename(z_path)
            copy(z_path, os.path.join(dst, z_file))
            decompress_file(os.path.join(dst, z_file))


def get_uri(dir_name: str) -> str:
    """
    Return the URI path for a directory.

    This allows files hosted on different file servers to have distinct locations.

    Adapted from atomate2.utils.path.get_uri.

    Parameters
    ----------
    dir_name : str or Path
        A directory name.

    Returns
    -------
    str
        Full URI path, e.g., "fileserver.host.com:/full/path/of/dir_name".
    """
    fullpath = Path(dir_name).absolute()
    hostname = socket.gethostname()
    with contextlib.suppress(socket.gaierror, socket.herror):
        hostname = socket.gethostbyaddr(hostname)[0]
    return f"{hostname}:{fullpath}"


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
