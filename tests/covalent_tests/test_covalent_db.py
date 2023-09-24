import os

import pytest
from maggma.stores import MemoryStore

from quacc import SETTINGS

try:
    import covalent as ct

    ct = ct if SETTINGS.WORKFLOW_ENGINE == "covalent" else None
except ImportError:
    ct = None


@pytest.mark.skipif(
    ct is None or os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions with Covalent",
)
def test_covalent_to_db():
    from quacc.utils.db import covalent_to_db

    store = MemoryStore(collection_name="db1")
    covalent_to_db(store)
    with store:
        count1 = store.count()
    assert count1 > 0

    with pytest.warns(UserWarning):
        covalent_to_db(store, dispatch_ids=["bad-value"])
    assert store.count() == count1

    store = MemoryStore(collection_name="db2")
    covalent_to_db(store, results_dir=ct.get_config()["dispatcher"]["results_dir"])
    with store:
        count2 = store.count()
    assert count2 == count1

    with pytest.raises(ValueError):
        covalent_to_db(store, dispatch_ids=["bad-value"], results_dir="bad-value")


@pytest.mark.skipif(
    ct is None or os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions with Covalent",
)
def test_covalent_db_tutorial():
    # Connect to the database

    store = MemoryStore()

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
