"""Utility functions for dealing with dictionaries."""
from __future__ import annotations
import copy
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any
def is_non_pickleable(obj):
    """
    Check if an object is non-pickleable.
    """
    try:
        copy.deepcopy(obj)
        return False
    except (TypeError, ValueError):
        return True
def handle_non_pickleable(obj):
    """
    Handle non-pickleable objects in a custom way.
    """
    if isinstance(obj, dict):
        # Handle non-pickleable dictionaries by creating a new one
        return {key: custom_deep_copy(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        # Handle non-pickleable lists by creating a new one
        return [custom_deep_copy(item) for item in obj]
    elif isinstance(obj, tuple):
        # Handle non-pickleable tuples by converting them to lists
        return [custom_deep_copy(item) for item in obj]
    elif isinstance(obj, set):
        # Handle non-pickleable sets by converting them to lists
        return [custom_deep_copy(item) for item in obj]
    elif isinstance(obj, CustomNonPickleableObject):
        # Handle custom non-pickleable objects as needed
        return obj.handle_special_case()
    else:
        # For other non-pickleable objects, return the original object
        return obj
    
class CustomNonPickleableObject:
    def __init__(self, data):
        self.data = data

    def handle_special_case(self):
        # Implement custom handling for this specific non-pickleable object
        return f"CustomNonPickleableObject({self.data})"

def custom_deep_copy(obj):
    if is_non_pickleable(obj):
        return handle_non_pickleable(obj)
    else:
        return copy.deepcopy(obj)

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
    merged = custom_deep_copy(dict1)  # Use custom_deep_copy here

    for key, value in dict2.items():
        if key in merged:
            if isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = merge_dicts(merged[key], value)
            else:
                merged[key] = custom_deep_copy(value)  # Use custom_deep_copy for values
        else:
            merged[key] = custom_deep_copy(value)  # Use custom_deep_copy for values

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
        old_dict = custom_deep_copy(merged)  # Use custom_deep_copy for copying the merged dictionary
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
