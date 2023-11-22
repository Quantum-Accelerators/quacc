"""Utility functions for dealing with dictionaries."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any


def merge_dicts(
    dict1: dict[str, Any] | None,
    dict2: dict[str, Any] | None,
    remove_nones: bool = True,
) -> dict[str, Any]:
    """
    Recursively merges two dictionaries. If one of the inputs is `None`, then it is
    treated as `{}`.

    This function should be used instead of the | operator when merging nested dictionaries,
    e.g. `{"a": {"b": 1}} | {"a": {"c": 2}}` will return `{"a": {"c": 2}}` whereas
    `merge_dicts({"a": {"b": 1}}, {"a": {"c": 2}})` will return `{"a": {"b": 1, "c": 2}}`.

    Parameters
    ----------
    dict1
        First dictionary
    dict2
        Second dictionary
    remove_nones
        If True, remove empty lists and dictionaries

    Returns
    -------
    dict
        Merged dictionary
    """
    dict1 = dict1 or {}
    dict2 = dict2 or {}
    merged = dict1.copy()

    for key, value in dict2.items():
        if key in merged:
            if isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = merge_dicts(merged[key], value)
            else:
                merged[key] = value
        else:
            merged[key] = value

    if remove_nones:
        merged = remove_dict_nones(merged)

    return merged


def merge_several_dicts(*args, remove_nones: bool = True) -> dict[str, Any]:
    """
    Recursively merge several dictionaries, taking the latter in the list as higher preference.

    Parameters
    ----------
    *args
        Dictionaries to merge
    remove_nones
        If True, remove empty lists and dictionaries

    Returns
    -------
    dict
        Merged dictionary
    """
    old_dict = args[0]
    for i in range(len(args) - 1):
        merged = merge_dicts(old_dict, args[i + 1], remove_nones=remove_nones)
        old_dict = merged.copy()
    return merged


def remove_dict_nones(start_dict: dict[str, Any]) -> dict[str, Any]:
    """
    For a given dictionary, recursively remove all items that are None.

    Parameters
    ----------
    start_dict
        Dictionary to clean

    Returns
    -------
    dict
        Cleaned dictionary
    """

    if isinstance(start_dict, dict):
        return {k: remove_dict_nones(v) for k, v in start_dict.items() if v is not None}
    return (
        [remove_dict_nones(v) for v in start_dict]
        if isinstance(start_dict, list)
        else start_dict
    )


def sort_dict(start_dict: dict[str, Any]) -> dict[str, Any]:
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
        k: sort_dict(v) if isinstance(v, dict) else v
        for k, v in sorted(start_dict.items())
    }
