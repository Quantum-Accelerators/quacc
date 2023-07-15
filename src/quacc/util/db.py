"""Utility functions for interfacing with databases"""
from __future__ import annotations

import os
import uuid
import warnings

import covalent as ct
from covalent._shared_files.exceptions import MissingLatticeRecordError
from maggma.core import Store
from monty.json import MontyDecoder


def covalent_to_db(
    store: Store, dispatch_ids: list[str] | None = None, results_dir: str | None = None
) -> None:
    """
    Store the results of a Covalent database in a user-specified Maggma Store

    Parameters
    ----------
    store
        The Maggma Store object to store the results in
    dispatch_ids
        Dispatch ID to store. If None, all dispatch IDs in the results_dir will be stored
    results_dir
        The Covalent results_dir to pull if dispatch_ID is None. If None, the results_dir from ct.get_config() will be used

    Returns
    -------
    None
    """

    if dispatch_ids and results_dir:
        raise ValueError("Cannot specify both dispatch_id and results_dir")
    dispatch_ids = dispatch_ids or []

    # Get the dispatch IDs
    if not dispatch_ids:
        if results_dir:
            dispatch_ids = os.listdir(results_dir)
        else:
            config_results_dir = ct.get_config()["dispatcher"]["results_dir"]
            dispatch_ids = os.listdir(config_results_dir)

    # Populate the docs
    docs = []
    for d_id in dispatch_ids:
        try:
            result_obj = ct.get_result(d_id)
        except MissingLatticeRecordError:
            warnings.warn(f"Could not find dispatch_id: {d_id}", UserWarning)
            continue
        if result_obj and result_obj.status == "COMPLETED":
            docs.append({"dispatch_id": d_id, "result": result_obj.result})

    # Store the results
    if docs:
        with store:
            store.update(docs, key="dispatch_id")


def results_to_db(store: Store | str, results: dict | list[dict]) -> None:
    """
    Store the results of a quacc recipe in a user-specified Maggma Store.
    A UUID will be generated for each entry.

    Parameters
    ----------
    store
        The Maggma Store object to store the results in or a str representation
        of a Maggma Store (taken from `.to_json()`)
    results
        The output summary dictionary or list of dictionaries from a quacc recipe

    Returns
    -------
    None
    """
    if isinstance(store, str):
        store = MontyDecoder().decode(store)

    if isinstance(results, dict):
        results = [results]
    for result in results:
        result["uuid"] = str(uuid.uuid4())

    with store:
        store.update(results, key="uuid")
