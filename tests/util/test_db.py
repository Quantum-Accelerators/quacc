import os

import covalent as ct
import pytest
from ase.build import bulk
from maggma.stores import MemoryStore, MontyStore

from quacc.recipes.emt.core import static_job
from quacc.util.db import covalent_to_db, results_to_db


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_covalent_to_db():
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


def test_results_to_db():
    atoms = bulk("Cu")
    output = static_job(atoms)
    store = MemoryStore(collection_name="db3")
    results_to_db(store, output)
    assert store.count() == 1


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_monty_db():
    store = MontyStore("quacc_collection", database_path=".", database_name="quacc_db")
    covalent_to_db(store, results_dir=ct.get_config()["dispatcher"]["results_dir"])
    count = store.count()
    assert count == 1

    atoms = bulk("Cu")
    output = static_job(atoms)
    results_to_db(store, output)
    assert store.count() == 2
