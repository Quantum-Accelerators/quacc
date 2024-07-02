"""Utility functions for file and path handling."""

from __future__ import annotations

import contextlib
import logging
import os
import socket
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from random import randint
from shutil import copy
from typing import TYPE_CHECKING

from monty.io import zopen
from monty.os.path import zpath
from monty.shutil import copy_r, decompress_dir, decompress_file
from ruamel.yaml import YAML

if TYPE_CHECKING:
    from typing import Any

    from quacc.types import Filenames, SourceDirectory


logger = logging.getLogger(__name__)


def check_logfile(logfile: str | Path, check_str: str) -> bool:
    """
    Check if a logfile has a given string (case-insensitive).
    The compression suffix, e.g. `.gz`, is automatically handled
    and does not need to be specified.

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
    logfile_path = Path(logfile).expanduser()
    zlog = Path(zpath(str(logfile_path)))
    with zopen(zlog, "r") as f:
        for line in f:
            clean_line = line if isinstance(line, str) else line.decode("utf-8")
            if check_str.lower() in clean_line.lower():
                return True
    return False


def copy_decompress_files(
    source_directory: SourceDirectory,
    filenames: Filenames,
    destination_directory: str | Path,
) -> None:
    """
    Copy and decompress `filenames` from the `source_directory` to the `destination`
    directory.

    For example, to copy the files `CHGCAR` and `WAVECAR` from the `source_directory` to
    the `destination` directory, use the following:

    ```python
    copy_decompress_files(
        source_directory="/path/to/source",
        filenames=["CHGCAR", "WAVECAR"],
        destination="/path/to/destination",
    )
    ```

    This function also supports glob patterns for any of the entries within `filenames`.

    For example, to copy and decompress all files in the `source_directory` with the
    extension `.gz` to the `destination` directory, use the following:

    ```python
    copy_decompress_files(
        source_directory="/path/to/source",
        filenames=["*.gz"],
        destination="/path/to/destination",
    )
    ```

    If a directory is specified in `filenames`, that directory and its contents will be
    copied and decompressed to the `destination` directory.

    For example, to recursively copy the entire directory `prior_run` and decompress its
    files from the `source_directory` to the `destination` directory, use the following:

    ```python
    copy_decompress_files(
        source_directory="/path/to/source",
        filenames=["prior_run"],
        destination="/path/to/destination",
    )
    ```

    Sometimes, you may want to copy a directory but only keep some of the files for the
    sake of saving space. In other words, you want to retain the tree structure of the
    files with respect to some parent directory. To do this, you can specify the
    tree to retain in the `filenames` argument. For example, to copy and decompress the
    files `prior_run/CHGCAR` and `prior_run/WAVECAR` from the `source_directory` to the
    `destination` directory while the tree structure, use the following:

    ```python
    copy_decompress_files(
        source_directory="/path/to/source",
        filenames=["prior_run/CHGCAR", "prior_run/WAVECAR"],
        destination="/path/to/destination",
    )
    ```

    Parameters
    ----------
    source_directory
        Directory to copy files from.
    filenames
        Files to copy and decompress. Glob patterns are supported.
    destination_directory
        Destination directory.

    Returns
    -------
    None
    """
    source_directory = Path(source_directory).expanduser()
    destination_directory = Path(destination_directory).expanduser()

    if not isinstance(filenames, list):
        filenames = [filenames]

    for f in filenames:
        globs_found = list(source_directory.glob(str(f)))
        if not globs_found:
            logger.warning(f"Cannot find file {f} in {source_directory}")
        for source_filepath in globs_found:
            destination_filepath = destination_directory / source_filepath.relative_to(
                source_directory
            )
            Path(destination_filepath.parent).mkdir(parents=True, exist_ok=True)

            if source_filepath.is_symlink():
                continue
            if source_filepath.is_file():
                copy(source_filepath, destination_filepath)
                decompress_file(destination_filepath)
            elif source_filepath.is_dir():
                copy_r(source_filepath, destination_filepath)
                decompress_dir(destination_filepath)


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
    config = YAML().load(yaml_path)

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


def find_recent_logfile(directory: Path | str, logfile_extensions: str | list[str]):
    """
    Find the most recent logfile in a given directory.

    Parameters
    ----------
    directory
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
    for f in Path(directory).expanduser().iterdir():
        f_path = Path(directory, f)
        for ext in logfile_extensions:
            if ext in str(f) and f_path.stat().st_mtime > mod_time:
                mod_time = f_path.stat().st_mtime
                logfile = f_path.resolve()
    return logfile


def get_uri(directory: str | Path) -> str:
    """
    Return the URI path for a directory.

    This allows files hosted on different file servers to have distinct
    locations.

    Adapted from Atomate2.

    Parameters
    ----------
    directory
        A directory name.

    Returns
    -------
    str
        Full URI path, e.g., "fileserver.host.com:/full/path/of/dir_name".
    """
    fullpath = Path(directory).expanduser().resolve()
    hostname = socket.gethostname()
    with contextlib.suppress(socket.gaierror, socket.herror):
        hostname = socket.gethostbyaddr(hostname)[0]
    return f"{hostname}:{fullpath}"


def safe_decompress_dir(path: str | Path) -> None:
    """
    Recursively decompresses all files in a directory.
    This is a wrapper around the `decompress_file` function.


    Args:
        path (str | Path): Path to parent directory.
    """
    path = Path(path)
    for parent, _, files in os.walk(path):
        for f in files:
            try:
                decompress_file(Path(parent, f))
            except FileNotFoundError:
                logger.debug(f"Cannot find {f} in {parent}. Skipping.")
