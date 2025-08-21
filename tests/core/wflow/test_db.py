from __future__ import annotations

import pytest

pytest.importorskip("maggma")
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.wflow_tools.db import results_to_db


def test_results_to_db():
    atoms = bulk("Cu")
    store = MemoryStore(collection_name="db3")
    results_to_db(store, {"atoms": atoms})
    with store:
        assert store.count() == 1
        assert store.query_one().get("uuid")
