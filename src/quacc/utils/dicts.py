"""
Utility functions for dealing with dictionaries
"""
from __future__ import annotations


def merge_dicts(
    dict1: dict | None, dict2: dict | None, remove_nones: bool = True
) -> dict:
    """
    Recursively merges two dictionaries. If one the inputs are `None`, then
    it is treated as `{}`.

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


def remove_dict_nones(start_dict: dict) -> dict:
    """
    For a given dictionary, recursively remove all items that are None

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


def sort_dict(start_dict: dict) -> dict:
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
