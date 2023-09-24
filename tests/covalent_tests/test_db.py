import os

import pytest
from maggma.stores import MemoryStore

try:
    import covalent as ct

except ImportError:
    ct = None


@pytest.mark.skipif(
    ct is None,
    reason="This test requires covalent to be the workflow engine",
)
def test_covalent_to_db():
    from quacc.utils.db import covalent_to_db

    store = MemoryStore(collection_name="db1")
    covalent_to_db(store)
    count1 = store.count()
    assert count1 > 0

    with pytest.warns(UserWarning):
        covalent_to_db(store, dispatch_ids=["bad-value"])
    assert store.count() == count1

    store = MemoryStore(collection_name="db2")
    covalent_to_db(store, results_dir=ct.get_config()["dispatcher"]["results_dir"])
    count2 = store.count()
    assert count2 == count1

    with pytest.raises(ValueError):
        covalent_to_db(store, dispatch_ids=["bad-value"], results_dir="bad-value")
