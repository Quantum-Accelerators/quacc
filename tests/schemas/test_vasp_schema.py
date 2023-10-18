import pytest


@pytest.fixture()
def run1():
    from pathlib import Path

    FILE_DIR = Path(__file__).resolve().parent

    return FILE_DIR / "vasp_run1"


def mock_bader_analysis(*args, **kwargs):
    return {
        "min_dist": [1.0] * 16,
        "charge": [1.0] * 16,
        "atomic_volume": [1.0] * 16,
        "vacuum_charge": 1.0,
        "vacuum_volume": 1.0,
        "bader_version": 1.0,
        "reference_used": [0.0] * 16,
        "charge_transfer": [-1.0] * 16,
        "magmom": [0.0] * 16,
    }


def test_vasp_summarize_run(run1):
    import os

    from ase.io import read
    from maggma.stores import MemoryStore
    from monty.json import MontyDecoder, jsanitize

    from quacc.calculators.vasp import Vasp
    from quacc.schemas.vasp import vasp_summarize_run

    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    calc = atoms.calc
    results = vasp_summarize_run(
        atoms,
        dir_path=run1,
    )
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["output"]["energy"] == -33.15807349

    # Make sure default dir works
    cwd = os.getcwd()
    os.chdir(run1)
    vasp_summarize_run(atoms)
    os.chdir(cwd)

    # Test DB
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    store = MemoryStore()
    vasp_summarize_run(atoms, dir_path=run1, store=store)
    assert store.count() == 1

    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = vasp_summarize_run(atoms, dir_path=run1, additional_fields={"test": "hi"})
    assert results.get("test") == "hi"

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0}
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = vasp_summarize_run(atoms, dir_path=run1)
    results_atoms = results["atoms"]
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results_atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0, "magmoms": [2.0] * len(atoms)}
    results = vasp_summarize_run(atoms, dir_path=run1)
    results_atoms = results["atoms"]

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results_atoms.get_initial_magnetic_moments().tolist() == [2.0] * len(atoms)
    assert results_atoms.calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, "CONTCAR.gz"))
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0}
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = vasp_summarize_run(atoms, dir_path=run1, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    results_atoms = results["atoms"]
    assert results_atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_bader_run(monkeypatch, run1, tmpdir):
    from shutil import copytree, move

    from ase.io import read

    from quacc.schemas.vasp import vasp_summarize_run

    monkeypatch.setattr(
        "quacc.schemas.vasp.bader_analysis_from_path",
        mock_bader_analysis,
    )
    tmpdir.chdir()

    p = tmpdir / "vasp_run"
    copytree(run1, p)

    move(p / "garbled_pot", p / "POTCAR")

    for f in ["CHGCAR.gz", "AECCAR0.gz", "AECCAR2.gz"]:
        with open(p / f, "w") as w:
            w.write("test")

    # Make sure Bader works
    atoms = read(str(p / "OUTCAR.gz"))
    results = vasp_summarize_run(atoms, dir_path=p, run_bader=True)
    struct = results["output"]["structure"]
    assert struct.site_properties["bader_charge"] == [1.0] * len(atoms)
    assert struct.site_properties["bader_spin"] == [0.0] * len(atoms)


def test_no_bader(run1, tmpdir):
    from ase.io import read

    from quacc.schemas.vasp import vasp_summarize_run

    tmpdir.chdir()

    atoms = read(run1 / "OUTCAR.gz")
    with pytest.warns(UserWarning):
        vasp_summarize_run(atoms, dir_path=run1, run_bader=True)
