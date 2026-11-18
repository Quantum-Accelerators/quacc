"""Utilities for working with lists."""

from __future__ import annotations


def merge_list_params(
    *lists: list[str] | None,
    removal_prefix: str | None = "#",
    case_insensitive: bool = True,
) -> list[str]:
    """
    Merge list string parameters, taking list2 as the priority. Also removes any entries
    that start with `removal_prefix` from the final list.

    Parameters
    ----------
    *lists
        Lists to merge, with the latter taking priority.
    removal_prefix
        Prefix to use to remove an entry from the final list.
    case_insensitive
        Whether to perform case-insensitive comparisons when identifying
        duplicate entries and entries to remove. The case of each entry is
        preserved in the returned list, since entries can contain
        case-sensitive content such as file paths.

    Returns
    -------
    list
        Merged list
    """
    merged: dict[str, str] = {}
    for list_ in lists:
        for item in list_ or []:
            key = item.lower() if case_insensitive else item
            merged[key] = item

    if removal_prefix:
        removal_keys = [key for key in merged if key.startswith(removal_prefix)]
        for key in removal_keys:
            merged.pop(key, None)
            merged.pop(key[len(removal_prefix) :], None)

    return [merged[key] for key in sorted(merged)]
