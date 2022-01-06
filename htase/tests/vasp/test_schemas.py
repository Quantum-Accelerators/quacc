import os
from htase.schemas.vasp.summarize import get_results
from htase.util.calc import cache_calc
from ase.io import read
from ase.io.jsonio import encode
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "run1")

# needs more tests
def test_summarize():
    atoms = read(os.path.join(run1, "POSCAR.gz"))
    results = get_results(dir_path=run1)

    results = get_results(atoms=atoms, dir_path=run1)
    assert results.get("atoms", None) is not None and results["atoms"] == encode(atoms)

    results = get_results(atoms=atoms, dir_path=run1)
    atoms = cache_calc(results["atoms"])
    assert (
        results["atoms"].get("results", None) is not None
        and results["atoms"]["results"].get("calc0", None) is not None
    )
