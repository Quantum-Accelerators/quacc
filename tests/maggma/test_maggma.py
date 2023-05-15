import os

import pytest
from maggma.stores import MemoryStore

try:
    import maggma as mg
except ImportError:
    mg = None
try:
    import covalent as ct
except ImportError:
    ct = None


@pytest.mark.skipif(
    mg is None or ct is None,
    reason="Maggma and Covalent must be installed",
)
def test_tutorial():
    # Connect to the database
    store = MemoryStore()
    store.connect()

    # Fetch the results
    results_dir = ct.get_config()["dispatcher"]["results_dir"]
    docs = []
    for i, dispatch_id in enumerate(os.listdir(results_dir)):
        result = ct.get_result(dispatch_id).result
        docs.append({"dispatch_id": dispatch_id, "result": result})
        if i > 2:
            break

    # Store the results
    with store:
        store.update(docs, key="dispatch_id")

    # Close the database connection
    store.close()
