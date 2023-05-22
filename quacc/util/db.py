import os
import warnings

import covalent as ct
from covalent._shared_files.exceptions import MissingLatticeRecordError
from maggma.core import Store


def covalent_to_db(
    store: Store = None, dispatch_id: str = None, results_dir: str = None
):
    """
    Store the results of a Covalent dispatch in a database

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
            result = ct.get_result(dispatch_id).result
        except MissingLatticeRecordError:
            warnings.warn(f"Could not find dispatch_id: {dispatch_id}", UserWarning)
        docs.append({"dispatch_id": dispatch_id, "result": result})

    # Store the results
    store.connect()
    with store:
        store.update(docs, key="dispatch_id")

    # Close the database connection
    store.close()
