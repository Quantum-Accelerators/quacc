"""Utility functions for dealing with dictionaries."""

from __future__ import annotations

from collections.abc import MutableMapping
from copy import deepcopy
from logging import getLogger
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

LOGGER = getLogger(__name__)


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
        old_dict = deepcopy(merged)
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
    merged = deepcopy(dict1)
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
