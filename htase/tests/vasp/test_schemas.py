import os
from htase.schemas.vasp.summarize import get_results
from htase.util.calc import cache_calc
from htase.calculators.vasp import SmartVasp
from ase.io import read
from ase.io.jsonio import encode, decode
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "run1")


def test_summarize():
    atoms = read(os.path.join(run1, "POSCAR.gz"))
    results = get_results(dir_path=run1)
    # assert something here

    results = get_results(atoms=atoms, dir_path=run1)
    assert results.get("atoms", None) is not None and results["atoms"] == encode(atoms)

    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0}
    results = get_results(atoms=atoms, dir_path=run1)
    atoms = decode(cache_calc(results["atoms"]))
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None)
    assert atoms.info["results"]["calc0"].get("energy", None) == -1.0
