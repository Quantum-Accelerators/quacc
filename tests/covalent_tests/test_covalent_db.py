import os

import pytest
from maggma.stores import MemoryStore

try:
    import covalent as ct

except ImportError:
    ct = None


@pytest.mark.skipif(
    ct is None or os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions with Covalent",
)
def test_covalent_to_db():
    import covalent as ct

    from quacc.utils.db import covalent_to_db

    @ct.electron
    def add(a, b):
        return a + b

    @ct.lattice
    def workflow(a, b):
        return add(a, b)

    dispatch_id = ct.dispatch(workflow)(1, 2)
    ct.get_result(dispatch_id, wait=True)
    dispatch_id = ct.dispatch(workflow)(1, 2)
    ct.get_result(dispatch_id, wait=True)

    store = MemoryStore(collection_name="db1")
    covalent_to_db(store)
    with store:
        count1 = store.count()
    assert count1 >= 2

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
