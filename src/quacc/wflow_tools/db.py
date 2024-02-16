"""Utility functions for interfacing with databases."""

from __future__ import annotations

import uuid
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

    from maggma.core import Store


def results_to_db(store: Store, results: dict[str, Any] | list[dict]) -> None:
    """
    Store the results of a quacc recipe in a user-specified Maggma Store. A UUID will be
    generated for each entry.

    Parameters
    ----------
    store
        The Maggma Store object to store the results in
    results
        The output summary dictionary or list of dictionaries from a quacc
        recipe

    Returns
    -------
    None
    """

    if not isinstance(results, list):
        results = [results]

    for result in results:
        result["uuid"] = str(uuid.uuid4())

    with store:
        store.update(results, key="uuid")
