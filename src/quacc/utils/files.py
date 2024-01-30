"""Utility functions for file and path handling."""
from __future__ import annotations

import contextlib
import socket
import warnings
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from random import randint
from shutil import copy
from typing import TYPE_CHECKING

import yaml
from monty.io import zopen
from monty.os.path import zpath
from monty.shutil import decompress_file

if TYPE_CHECKING:
    from typing import Any


def check_logfile(logfile: str, check_str: str) -> bool:
    """
    Check if a logfile has a given string (case-insensitive).

    Parameters
    ----------
    logfile
        Path to the logfile.
    check_str
        String to check for.

    Returns
    -------
    bool
        True if the string is found in the logfile, False otherwise.
    """
    zlog = Path(zpath(logfile)).expanduser()
    with zopen(zlog, "r") as f:
        for line in f:
            clean_line = line if isinstance(line, str) else line.decode("utf-8")
            if check_str.lower() in clean_line.lower():
                return True
    return False


def copy_decompress_files(
    source_files: list[str | Path], destination: str | Path
) -> None:
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
        f_path = Path(zpath(f)).expanduser()
        if f_path.is_symlink():
            continue
        if f_path.is_file():
            copy(f_path, Path(destination, f_path.name))
            decompress_file(Path(destination, f_path.name))
        elif f_path.is_dir():
            copy_decompress_files_from_dir(f_path, destination)
        else:
            warnings.warn(f"Cannot find file {f_path}", UserWarning)


def copy_decompress_tree(
    source_files: dict[str, str | Path | list[str | Path]], destination: str | Path
) -> None:
    """
    Copy and decompress files from source to destination. This function respects the
    directory tree.

    Parameters
    ----------
    source_files
        Dict, key is the base_dir, values are the tree to
        be respected
    destination
        Destination directory.

    Returns
    -------
    None
    """

    # Work with glob pattern, work if the glob pattern return nothing
    for _base, _tree in source_files.items():
        base = Path(_base).expanduser()

        abs_files, rel_files = [], []

        tree = [_tree] if not isinstance(_tree, list) else _tree

        for f in tree:
            glob_found = list(base.glob(f))
            abs_files.extend(glob_found)
            rel_files.extend([i.relative_to(base) for i in glob_found])

        for abs_f, rel_f in zip(abs_files, rel_files):
            Path(destination, rel_f.parent).mkdir(parents=True, exist_ok=True)
            copy_decompress_files([abs_f], Path(destination, rel_f.parent))


def copy_decompress_files_from_dir(source: str | Path, destination: str | Path) -> None:
    """
    Copy and decompress files recursively from source to destination.

    Parameters
    ----------
    source
        Directory to walk and copy files from.
    destination
        Destination directory.

    Returns
    -------
    None
    """
    src, dst = Path(source).expanduser(), Path(destination).expanduser()

    if src.is_dir():
        for f in src.iterdir():
            if f.resolve() == dst.resolve() or f.is_symlink():
                continue
            if f.is_file():
                copy_decompress_files([f], dst)
            elif f.is_dir():
                (dst / f.name).mkdir(exist_ok=True)
                copy_decompress_files_from_dir(src / f, dst / f.name)
    else:
        warnings.warn(f"Cannot find {src}", UserWarning)


def make_unique_dir(
    base_path: Path | str | None = None, prefix: str | None = None
) -> Path:
    """
    Make a directory with a unique name.

    Parameters
    ----------
    base_path
        Path to the base directory.
    prefix
        Prefix to add to the directory name.

    Returns
    -------
    Path
        Path to the job directory.
    """
    time_now = datetime.now(timezone.utc).strftime("%Y-%m-%d-%H-%M-%S-%f")
    if prefix is None:
        prefix = ""
    job_dir = Path(f"{prefix}{time_now}-{randint(10000, 99999)}")
    if base_path:
        job_dir = Path(base_path, job_dir)
    job_dir.mkdir(parents=True)

    return job_dir


def load_yaml_calc(yaml_path: str | Path) -> dict[str, Any]:
    """
    Loads a YAML file containing calculator settings. This YAML loader looks for a
    special flag "parent" in the YAML file. If this flag is present, the YAML file
    specified in the "parent" flag is loaded and its contents are inherited by the child
    YAML file.

    Parameters
    ----------
    yaml_path
        Path to the YAML file.

    Returns
    -------
    dict
        The calculator configuration (i.e. settings).
    """

    yaml_path = Path(yaml_path).expanduser()

    if yaml_path.suffix != ".yaml":
        yaml_path = yaml_path.with_suffix(f"{yaml_path.suffix}.yaml")

    if not yaml_path.exists():
        msg = f"Cannot find {yaml_path}"
        raise FileNotFoundError(msg)

    # Load YAML file
    with yaml_path.open() as stream:
        config = yaml.safe_load(stream)

    # Inherit arguments from any parent YAML files but do not overwrite those in
    # the child file.
    for config_arg in deepcopy(config):
        if "parent" in config_arg.lower():
            yaml_parent_path = Path(yaml_path).parent / Path(config[config_arg])
            parent_config = load_yaml_calc(yaml_parent_path)

            for k, v in parent_config.items():
                if k not in config:
                    config[k] = v
                else:
                    v_new = parent_config.get(k, {})
                    for kk, vv in v_new.items():
                        if kk not in config[k]:
                            config[k][kk] = vv

            del config[config_arg]

    return config


def find_recent_logfile(dir_name: Path | str, logfile_extensions: str | list[str]):
    """
    Find the most recent logfile in a given directory.

    Parameters
    ----------
    dir_name
        The path to the directory to search
    logfile_extensions
        The extension (or list of possible extensions) of the logfile to search
        for. For an exact match only, put in the full file name.

    Returns
    -------
    logfile
        The path to the most recent logfile with the desired extension
    """
    mod_time = 0.0
    logfile = None
    if isinstance(logfile_extensions, str):
        logfile_extensions = [logfile_extensions]
    for f in Path(dir_name).expanduser().iterdir():
        f_path = Path(dir_name, f)
        for ext in logfile_extensions:
            if ext in str(f) and f_path.stat().st_mtime > mod_time:
                mod_time = f_path.stat().st_mtime
                logfile = f_path.resolve()
    return logfile


def get_uri(dir_name: str | Path) -> str:
    """
    Return the URI path for a directory.

    This allows files hosted on different file servers to have distinct
    locations.

    Adapted from Atomate2.

    Parameters
    ----------
    dir_name
        A directory name.

    Returns
    -------
    str
        Full URI path, e.g., "fileserver.host.com:/full/path/of/dir_name".
    """
    fullpath = Path(dir_name).expanduser().resolve()
    hostname = socket.gethostname()
    with contextlib.suppress(socket.gaierror, socket.herror):
        hostname = socket.gethostbyaddr(hostname)[0]
    return f"{hostname}:{fullpath}"
