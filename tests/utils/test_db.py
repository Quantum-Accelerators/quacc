from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import static_job
from quacc.utils.db import results_to_db


def test_results_to_db():
    atoms = bulk("Cu")
    output = static_job(atoms)
    store = MemoryStore(collection_name="db3")
    results_to_db(store, output)
    assert store.count() == 1
