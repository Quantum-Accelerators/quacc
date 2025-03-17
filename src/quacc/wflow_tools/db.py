"""Utility functions for interfacing with databases."""

from __future__ import annotations

import uuid
from functools import wraps
from typing import TYPE_CHECKING, Any

from monty.json import jsanitize

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from maggma.core import Store


def store_wrapper(func: Callable[..., Any], store: Store | None) -> Callable[..., Any]:
    """
    Wrap a function to store the results in the database.

    Parameters
    ----------
    func
        The function to wrap.
    settings
        The Quacc settings.

    Returns
    -------
    Callable[..., Any]
        The wrapped function.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)

        if isinstance(result, dict) and store:
            results_to_db(store, result)

        return result

    return wrapper


def results_to_db(store: Store, result: dict[str, Any]) -> None:
    """
    Store the results of a quacc recipe in a user-specified Maggma Store. A UUID will be
    generated for each entry.

    Parameters
    ----------
    store
        The Maggma Store object to store the results in
    results
        The output summary dictionary from a quacc recipe

    Returns
    -------
    None
    """

    sanitized_result = jsanitize(result, enum_values=True, recursive_msonable=True)
    sanitized_result["uuid"] = str(uuid.uuid4())

    with store:
        store.update(sanitized_result, key="uuid")
