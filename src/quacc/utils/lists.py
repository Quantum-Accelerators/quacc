"""Utilities for working with lists."""

from __future__ import annotations


def merge_list_params(
    *lists: list[str], removal_prefix: str | None = "#", case_insensitive: bool = True
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
        Whether to perform case-insensitive comparisons

    Returns
    -------
    list
        Merged list
    """
    merged_list = []
    lists_ = [list_ for list_ in lists if list_]
    for list_ in lists_:
        for item in list_:
            item_ = item
            if case_insensitive:
                item_ = item.lower()
            if item_ not in merged_list:
                merged_list.append(item_)

    if removal_prefix:
        items_to_remove1 = [
            item[1:] for item in merged_list if item.startswith(removal_prefix)
        ]
        items_to_remove2 = [
            item for item in merged_list if item.startswith(removal_prefix)
        ]
        items_to_remove = items_to_remove1 + items_to_remove2
        for item in items_to_remove:
            if item in merged_list:
                merged_list.remove(item)

    merged_list.sort()
    return merged_list
