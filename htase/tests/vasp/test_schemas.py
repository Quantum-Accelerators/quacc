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
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    results = get_results(dir_path=run1)
    assert results["nsites"] == len(atoms)

    results = get_results(atoms=atoms, dir_path=run1)
    assert results.get("atoms", None) is not None and results["atoms"] == encode(atoms)

    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.info["test"] = "hi"
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar", "cow": 1}
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0, "magmoms": [2.0] * len(atoms)}
    results = get_results(atoms=atoms, dir_path=run1)
    atoms = decode(cache_calc(results["atoms"]))
    assert atoms.info.get("results", None) is not None
    assert atoms.info["results"].get("calc0", None) is not None
    assert atoms.info["results"]["calc0"].get("energy", None) == -1.0
    assert atoms.info["results"]["calc0"].get("magmoms", None) == [2.0] * len(atoms)
    assert atoms.info["results"]["calc0"].get("rundir", None) is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [2.0] * len(atoms)
    assert atoms.info["test"] == "hi"
    assert atoms.info["test_dict"] == {"hi": "there", "foo": "bar", "cow": 1}
    assert results.get("atoms_info") is not None
    assert results["atoms_info"].get("test", None) == "hi"
    assert results["atoms_info"]["test_dict"] == {"hi": "there", "foo": "bar", "cow": 1}
