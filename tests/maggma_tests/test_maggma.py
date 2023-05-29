import os

import covalent as ct
import pytest

try:
    import maggma
except ImportError:
    maggma = None


@pytest.mark.skipif(
    maggma is None,
    reason="This test requiers maggma",
)
def test_tutorial():
    from maggma.stores import MemoryStore

    # Connect to the database

    store = MemoryStore()
    store.connect()

    # Fetch the results
    results_dir = ct.get_config()["dispatcher"]["results_dir"]
    docs = []
    for i, dispatch_id in enumerate(os.listdir(results_dir)):
        result = ct.get_result(dispatch_id).result
        docs.append({"dispatch_id": dispatch_id, "result": result})
        if i >= 2:
            break

    # Store the results
    with store:
        store.update(docs, key="dispatch_id")

    # Close the database connection
    store.close()
