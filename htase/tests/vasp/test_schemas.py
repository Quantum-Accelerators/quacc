import os
from htase.schemas.vasp.summarize import get_results
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

    atoms.info = {"test": "hi", "test2": ["hi", "bye"]}
    results = get_results(atoms=atoms, dir_path=run1)
    assert results.get("test", None) == "hi"
    assert results.get("test2", None) == ["hi", "bye"]
