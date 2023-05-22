import covalent as ct
import pytest
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import static_job
from quacc.util.db import covalent_to_db, results_to_db


def test_covalent_to_db():
    store = MemoryStore(collection_name="db1")
    covalent_to_db(store)
    count1 = store.count()
    assert count1 > 0

    with pytest.warns(UserWarning):
        covalent_to_db(store, dispatch_id="bad-value")
    assert store.count() == count1

    store = MemoryStore(collection_name="db2")
    covalent_to_db(store, results_dir=ct.get_config()["dispatcher"]["results_dir"])
    count2 = store.count()
    assert count2 == count1

    with pytest.raises(ValueError):
        covalent_to_db(store, dispatch_id="bad-value", results_dir="bad-value")


def test_results_to_db():
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += [0.1, 0.1, 0.1]

    output = static_job(atoms)
    store = MemoryStore(collection_name="db3")
    results_to_db(store, output)
    assert store.count() == 1
