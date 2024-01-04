"""
Utilities for working with lists.
"""
from __future__ import annotations


def merge_list_params(
    *lists: list[str], removal_prefix: str = "#", case_sensitive: bool = False
) -> list[str]:
    """
    Merge list string parameters, taking list2 as the priority.
    Also removes any entries that start with `removal_prefix` from the final list.

    Parameters
    ----------
    *lists
        Lists to merge, with the latter taking priority.
    removal_prefix
        Prefix to use to remove an entry from the final list.
    case_sensitive
        Whether to perform case-sensitive comparisons

    Returns
    -------
    list
        Merged list
    """

    lists = [list_ or [] for list_ in lists]
    old_list = lists[0]
    for i in range(len(lists) - 1):
        for j, element in enumerate(lists[i + 1]):
            if element.startswith(removal_prefix):
                old_list = [item for item in old_list if item != element[1:]]
                lists[i + 1][j] = element[1:]
        merged_list = _merge_list_pair(
            old_list, lists[i + 1], case_sensitive=case_sensitive
        )
        old_list = merged_list

    merged_list.sort()
    return merged_list


def _merge_list_pair(
    list1: list[str] | None, list2: list[str] | None, case_sensitive: bool = False
) -> list[str]:
    """
    Merges two lists of strings, with the latter taking priority.

    Parameters
    ----------
    list1
        First list
    list2
        Second list
    case_sensitive
        Whether to perform case-sensitive comparisons

    Returns
    -------
    list
        Merged list
    """
    list1 = list1 or []
    list2 = list2 or []

    if not case_sensitive:
        list1 = [element.lower() for element in list1]
        list2 = [element.lower() for element in list2]

    return [element for element in list1 if element not in list2] + list2
