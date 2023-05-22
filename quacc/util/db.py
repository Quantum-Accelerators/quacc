from __future__ import annotations

import os
import uuid
import warnings

import covalent as ct
from covalent._shared_files.exceptions import MissingLatticeRecordError
from maggma.core import Store


def covalent_to_db(store: Store, dispatch_id: str = None, results_dir: str = None):
    """
    Store the results of a Covalent database in a user-specified Maggma Store

    Parameters
    ----------
    store
        The Maggma Store object to store the results in
    dispatch_id
        A specific dispatch ID to store. If None, all dispatch IDs in the results_dir will be stored
    results_dir
        The Covalent results_dir to pull if dispatch_ID is None. If None, the results_dir from ct.get_config() will be used

    Returns
    -------
    None
    """

    if dispatch_id and results_dir:
        raise ValueError("Cannot specify both dispatch_id and results_dir")

    # Get the dispatch IDs
    dispatch_ids = []
    if dispatch_id:
        dispatch_ids = [dispatch_id]
    else:
        if results_dir:
            dispatch_ids = os.listdir(results_dir)
        else:
            config_results_dir = ct.get_config()["dispatcher"]["results_dir"]
            dispatch_ids = os.listdir(config_results_dir)

    # Populate the docs
    docs = []
    for dispatch_id in dispatch_ids:
        try:
            result_obj = ct.get_result(dispatch_id)
        except MissingLatticeRecordError:
            warnings.warn(f"Could not find dispatch_id: {dispatch_id}", UserWarning)
            continue
        if result_obj.status == "COMPLETED":
            docs.append({"dispatch_id": dispatch_id, "result": result_obj.result})

    # Store the results
    if docs:
        store.connect()
        with store:
            store.update(docs, key="dispatch_id")
        print(f"Stored {len(docs)} results in your database.")
        store.close()


def results_to_db(results: dict | list[dict], store: Store):
    """
    Store the results of a Quacc recipe in a user-specified Maggma Store.
    A UUID will be generated for each entry.

    Parameters
    ----------
    results
        The output summary dictionary or list of dictionaries from a Quacc recipe
    store
        The Maggma Store object to store the results in

    Returns
    -------
    None
    """

    if isinstance(results, dict):
        results = [results]
    for result in results:
        result["uuid"] = str(uuid.uuid4())

    store.connect()
    with store:
        store.update(results, key="uuid")
    print(f"Stored {len(results)} results in your database.")
    store.close()
