"""
Utility functions for dealing with dictionaries
"""
from __future__ import annotations


def clean_dict(d: dict, remove_empties=False) -> dict:
    """
    For a given dictionary, recursively remove all items that are None
    or are empty lists/dicts, and then sort all entries alphabetically by key.

    Parameters
    ----------
    d
        Dictionary to clean
    remove_empties
        If True, remove empty lists and dictionaries

    Returns
    -------
    dict
        Cleaned and sorted dictionary
    """

    if remove_empties:
        d = remove_dict_empties(d)
    return sort_dict(d)


def remove_dict_empties(d: dict) -> dict:
    """
    For a given dictionary, recursively remove all items that are None
    or are empty lists/dicts.

    Parameters
    ----------
    d
        Dictionary to clean

    Returns
    -------
    dict
        Cleaned dictionary
    """

    if isinstance(d, dict):
        return {
            k: remove_dict_empties(v)
            for k, v in d.items()
            if v is not None and (not isinstance(v, (dict, list)) or len(v) != 0)
        }
    return [remove_dict_empties(v) for v in d] if isinstance(d, list) else d


def sort_dict(d: dict) -> dict:
    """
    For a given dictionary, recursively sort all entries alphabetically by key.

    Parameters
    ----------
    d
        Dictionary to sort

    Returns
    -------
    dict
        Sorted dictionary
    """

    return {k: sort_dict(v) if isinstance(v, dict) else v for k, v in sorted(d.items())}
