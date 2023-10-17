from pathlib import Path

import pytest

FILE_DIR = Path(__file__).resolve().parent

run1 = FILE_DIR / "gaussian_run1"
log1 = run1 / "Gaussian.log"


@pytest.fixture()
def cclib_obj():
    from cclib.io import ccread

    return ccread(log1)


def setup_module():
    import gzip
    import shutil

    p = FILE_DIR / "cclib_data"

    with gzip.open(p / "psi_test.cube.gz", "r") as f_in, open(
        p / "psi_test.cube", "wb"
    ) as f_out:
        shutil.copyfileobj(f_in, f_out)


def teardown_module():
    import os

    p = FILE_DIR / "cclib_data"

    if os.path.exists(p / "psi_test.cube"):
        os.remove(p / "psi_test.cube")


def bad_mock_cclib_calculate(*args, **kwargs):
    msg = "This is a bad run."
    raise ValueError(msg)


def test_cclib_summarize_run():
    import os

    from ase.io import read
    from maggma.stores import MemoryStore
    from monty.json import MontyDecoder, jsanitize

    from quacc.schemas.cclib import cclib_summarize_run

    # Make sure metadata is made
    atoms = read(log1)
    results = cclib_summarize_run(atoms, ".log", dir_path=run1)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["spin_multiplicity"] == 1
    assert results["natoms"] == 6
    assert results["results"]["metadata"].get("success", None) is True
    assert results["results"].get("energy", None) == pytest.approx(-5516.118738093933)
    assert "pull_request" in results["builder_meta"]

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)

    # Make sure default dir works
    cwd = os.getcwd()
    os.chdir(run1)
    cclib_summarize_run(atoms, ".log")
    os.chdir(cwd)

    # Test DB
    atoms = read(log1)
    store = MemoryStore()
    cclib_summarize_run(atoms, ".log", dir_path=run1, store=store)
    assert store.count() == 1

    # Make sure info tags are handled appropriately
    atoms = read(log1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = cclib_summarize_run(atoms, ".log", dir_path=run1)
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = cclib_summarize_run(atoms, ".log", dir_path=run1)

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results["atoms"].get_initial_magnetic_moments().tolist() == [2.0] * len(
        atoms
    )
    assert results["atoms"].calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, log1))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = cclib_summarize_run(atoms, ".log", dir_path=run1, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    assert results["atoms"].get_initial_magnetic_moments().tolist() == [3.14] * len(
        atoms
    )


def test_errors():
    from ase.build import bulk

    from quacc.calculators.vasp import Vasp
    from quacc.schemas.cclib import cclib_summarize_run

    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        cclib_summarize_run(atoms, ".log", dir_path=run1)

    calc = Vasp(atoms)
    atoms.calc = calc
    with pytest.raises(ValueError):
        cclib_summarize_run(atoms, ".log", dir_path=run1)


def test_cclib_taskdoc(tmpdir):
    import os

    from monty.json import MontyDecoder, jsanitize

    from quacc.schemas.cclib import _cclibTaskDocument

    tmpdir.chdir()

    p = FILE_DIR / "cclib_data"

    # Plain parsing of task doc. We do not check all cclib entries
    # because they will evolve over time. We only check the ones we have
    # added and some important ones.
    doc = _cclibTaskDocument.from_logfile(p, ".log.gz").dict()
    assert doc["additional_attributes"]["energy"] == pytest.approx(-4091.763)
    assert doc["natoms"] == 2
    assert doc["charge"] == 0
    assert doc["spin_multiplicity"] == 3
    assert doc["nelectrons"] == 16
    assert "schemas" in doc["dir_name"]
    assert "gau_testopt.log.gz" in doc["logfile"]
    assert doc.get("attributes") is not None
    assert doc["attributes"]["metadata"]["success"] is True
    assert doc["molecule_initial"][0].coords == pytest.approx([0, 0, 0])
    assert doc["molecule"][0].coords == pytest.approx([0.397382, 0.0, 0.0])
    assert doc["additional_attributes"]["homo_energies"] == pytest.approx(
        [-7.054007346511501, -11.618445074798501]
    )
    assert doc["additional_attributes"]["lumo_energies"] == pytest.approx(
        [4.2384453353880005, -3.9423854660440005]
    )
    assert doc["additional_attributes"]["homo_lumo_gaps"] == pytest.approx(
        [11.292452681899501, 7.6760596087545006]
    )
    assert doc["additional_attributes"]["min_homo_lumo_gap"] == pytest.approx(
        7.6760596087545006
    )

    # Now we will try two possible extensions, but we will make sure that
    # it fails because the newest log file (.txt) is not valid
    with open(p / "test.txt", "w") as f:
        f.write("I am a dummy log file")
    with pytest.raises(Exception) as e:
        doc = _cclibTaskDocument.from_logfile(p, [".log", ".txt"]).dict()
    os.remove(p / "test.txt")
    assert "Could not parse" in str(e.value)

    # Test a population analysis
    doc = _cclibTaskDocument.from_logfile(p, "psi_test.out", analysis="MBO").dict()
    assert doc["attributes"]["mbo"] is not None

    # Let's try with two analysis (also check case-insensitivity)
    doc = _cclibTaskDocument.from_logfile(
        p, "psi_test.out", analysis=["mbo", "density"]
    ).dict()
    assert doc["attributes"]["mbo"] is not None
    assert doc["attributes"]["density"] is not None

    # Test a population analysis that will fail
    doc = _cclibTaskDocument.from_logfile(p, ".log", analysis="MBO").dict()
    assert doc["attributes"]["mbo"] is None

    doc = _cclibTaskDocument.from_logfile(p, "psi_test.out", analysis=["Bader"]).dict()
    assert doc["attributes"]["bader"] is not None

    # Make sure storing the trajectory works
    doc = _cclibTaskDocument.from_logfile(p, ".log", store_trajectory=True).dict()
    assert len(doc["trajectory"]) == 7
    assert doc["trajectory"][0] == doc["molecule_initial"]
    assert doc["trajectory"][-1] == doc["molecule"]

    # Make sure additional fields can be stored
    doc = _cclibTaskDocument.from_logfile(
        p, ".log", additional_fields={"test": "hi"}
    ).dict()
    assert doc["test"] == "hi"

    with pytest.raises(FileNotFoundError):
        _cclibTaskDocument.from_logfile(p, "does_not_exists.txt")

    # test document can be jsanitized
    d = jsanitize(doc, enum_values=True)

    # and decoded
    MontyDecoder().process_decoded(d)


def test_cclib_calculate(tmpdir, cclib_obj):
    from quacc.schemas.cclib import _cclib_calculate

    tmpdir.chdir()

    with pytest.raises(ValueError):
        _cclib_calculate(cclib_obj, method="fakemethod")

    with pytest.raises(ValueError):
        _cclib_calculate(cclib_obj, method="bader")

    with pytest.raises(FileNotFoundError):
        _cclib_calculate(cclib_obj, method="bader", cube_file="does_not_exists.txt")

    with pytest.raises(FileNotFoundError):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "cclib_data" / "psi_test.cube",
            proatom_dir="does_not_exists",
        )

    with pytest.raises(ValueError):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "cclib_data" / "psi_test.cube",
        )

    with pytest.raises(Exception):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "cclib_data" / "psi_test.cube",
            proatom_dir=FILE_DIR / "cclib_data" / "psi_test.cube",
        )


def test_monkeypatches(tmpdir, monkeypatch, cclib_obj):
    from quacc.schemas.cclib import _cclib_calculate

    tmpdir.chdir()
    monkeypatch.setenv("PROATOM_DIR", str(FILE_DIR / "cclib_data" / "proatomdata"))
    with pytest.raises(FileNotFoundError):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "cclib_data" / "psi_test.cube",
        )

    monkeypatch.setattr("cclib.method.Bader.calculate", bad_mock_cclib_calculate)
    with pytest.warns(UserWarning):
        assert (
            _cclib_calculate(
                cclib_obj,
                method="bader",
                cube_file=FILE_DIR / "cclib_data" / "psi_test.cube",
            )
            is None
        )
