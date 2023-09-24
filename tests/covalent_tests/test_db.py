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
def test_results_to_db():
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

    dispatch_id2 = ct.dispatch(workflow)(1, 2)
    ct.get_result(dispatch_id2, wait=True)
    store = MemoryStore(collection_name="db")
    covalent_to_db(store, dispatch_ids=[dispatch_id, dispatch_id2])
    assert store.count() == 2
