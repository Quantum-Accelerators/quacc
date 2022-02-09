from typing import Dict, Any

# Functions in here should be ported to monty if practical


def merge_dicts(
    d1: Dict[str, Any],
    d2: Dict[str, Any],
    remove_none: bool = False,
    remove_false: bool = False,
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
    """
    d1_clean = {k.lower(): v for k, v in d1.items()}
    d2_clean = {k.lower(): v for k, v in d2.items()}
    d_merged = {**d1_clean, **d2_clean}
    if remove_none:
        d_merged = {k: v for k, v in d_merged.items() if v is not None}
    if remove_false:
        d_merged = {k: v for k, v in d_merged.items() if v is not False}
    return d_merged
