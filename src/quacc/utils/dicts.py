"""Utility functions for dealing with dictionaries."""

from __future__ import annotations

import logging
from collections.abc import MutableMapping
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

from monty.json import jsanitize
from monty.serialization import dumpfn

from quacc.wflow_tools.db import results_to_db

if TYPE_CHECKING:
    from typing import Any

    from maggma.stores import Store

LOGGER = logging.getLogger(__name__)


class Remove:
    """
    A sentinel class used in quacc to mark a key in a dictionary for removal.

    Note: This is more robust than using `None` as the sentinel value because
    `None` is a valid value for many keyword arguments.
    """

    def __init__(self):
        raise NotImplementedError(
            "Remove is a sentinel class and should not be instantiated."
        )


def recursive_dict_merge(
    *dicts: MutableMapping[str, Any] | None,
    remove_trigger: Any = Remove,
    verbose: bool = False,
) -> MutableMapping[str, Any]:
    """
    Recursively merge several dictionaries, taking the latter in the list as higher
    preference. Also removes any entries that have a value of `remove_trigger` from the
    final dictionary. If a `None` is provided, it is assumed to be `{}`.

    This function should be used instead of the | operator when merging nested dictionaries,
    e.g. `{"a": {"b": 1}} | {"a": {"c": 2}}` will return `{"a": {"c": 2}}` whereas
    `recursive_dict_merge({"a": {"b": 1}}, {"a": {"c": 2}})` will return `{"a": {"b": 1, "c": 2}}`.

    Parameters
    ----------
    *dicts
        Dictionaries to merge
    remove_trigger
        Value to that triggers removal of the entry
    verbose
        Whether to log warnings when overwriting keys

    Returns
    -------
    MutableMapping[str, Any]
        Merged dictionary
    """
    old_dict = dicts[0]
    for i in range(len(dicts) - 1):
        merged = _recursive_dict_pair_merge(old_dict, dicts[i + 1], verbose=verbose)
        old_dict = safe_dict_copy(merged)
    return remove_dict_entries(merged, remove_trigger=remove_trigger)


def _recursive_dict_pair_merge(
    dict1: MutableMapping[str, Any] | None,
    dict2: MutableMapping[str, Any] | None,
    verbose: bool = False,
) -> MutableMapping[str, Any]:
    """
    Recursively merges two dictionaries. If a `None` is provided, it is assumed to be `{}`.

    Parameters
    ----------
    dict1
        First dictionary
    dict2
        Second dictionary
    verbose
        Whether to log warnings when overwriting keys

    Returns
    -------
    dict
        Merged dictionary
    """
    dict1 = dict1 or ({} if dict1 is None else dict1.__class__())
    dict2 = dict2 or ({} if dict2 is None else dict2.__class__())
    merged = safe_dict_copy(dict1)

    for key, value in dict2.items():
        if key in merged:
            if isinstance(merged[key], MutableMapping) and isinstance(
                value, MutableMapping
            ):
                merged[key] = _recursive_dict_pair_merge(
                    merged[key], value, verbose=verbose
                )
            else:
                merged[key] = value
                if verbose:
                    LOGGER.warning(f"Overwriting key '{key}' to: '{merged[key]}'")
        else:
            merged[key] = value

    return merged


def safe_dict_copy(d: MutableMapping[str, Any] | dict[str, Any]) -> dict:
    """
    Safely copy a dictionary to account for deepcopy errors.

    Parameters
    ----------
    d
        Dictionary to copy

    Returns
    -------
    dict
        Copied dictionary
    """

    try:
        return deepcopy(dict(d))
    except Exception:
        return dict(d).copy()


def remove_dict_entries(
    start_dict: MutableMapping[str, Any], remove_trigger: Any
) -> MutableMapping[str, Any]:
    """
    For a given dictionary, recursively remove all items that are the `remove_trigger`.

    Parameters
    ----------
    start_dict
        Dictionary to clean
    remove_trigger
        Value to that triggers removal of the entry

    Returns
    -------
    dict
        Cleaned dictionary
    """
    if isinstance(start_dict, MutableMapping):
        return {
            k: remove_dict_entries(v, remove_trigger)
            for k, v in start_dict.items()
            if v is not remove_trigger
        }
    return (
        [remove_dict_entries(v, remove_trigger) for v in start_dict]
        if isinstance(start_dict, list)
        else start_dict
    )


def sort_dict(start_dict: MutableMapping[str, Any]) -> MutableMapping[str, Any]:
    """
    For a given dictionary, recursively sort all entries alphabetically by key.

    Parameters
    ----------
    start_dict
        Dictionary to sort

    Returns
    -------
    dict
        Sorted dictionary
    """
    return {
        k: sort_dict(dict(v)) if isinstance(v, MutableMapping) else v
        for k, v in sorted(start_dict.items())
    }


def clean_dict(start_dict: MutableMapping[str, Any]) -> MutableMapping[str, Any]:
    """
    Clean up a task document dictionary by removing all entries that are None and
    sorting the dictionary alphabetically by key.

    Parameters
    ----------
    start_dict
        Dictionary to clean

    Returns
    -------
    dict
        Cleaned dictionary
    """
    return sort_dict(remove_dict_entries(start_dict, None))


def finalize_dict(
    task_doc: dict,
    directory: str | Path | None = None,
    gzip_file: bool = True,
    store: Store | None = None,
) -> MutableMapping[str, Any]:
    """
    Finalize a schema by cleaning it and storing it in a database and/or file.

    Parameters
    ----------
    task_doc
        Dictionary representation of the task document.
    directory
        Directory where the results file is stored.
    gzip_file
        Whether to gzip the results file.
    store
        Maggma Store object to store the results in.

    Returns
    -------
    dict
        Cleaned task document
    """

    cleaned_task_doc = clean_dict(task_doc)
    if directory:
        if "tmp-quacc" in str(directory):
            raise ValueError("The directory should not be a temporary directory.")

        sanitized_schema = jsanitize(
            cleaned_task_doc, enum_values=True, recursive_msonable=True
        )
        dumpfn(
            sanitized_schema,
            Path(
                directory,
                "quacc_results.json.gz" if gzip_file else "quacc_results.json",
            ),
            fmt="json",
            indent=4,
        )

    if store:
        results_to_db(store, task_doc)

    return cleaned_task_doc
