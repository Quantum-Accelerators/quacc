
from __future__ import annotations

from typing import Any, Dict


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
