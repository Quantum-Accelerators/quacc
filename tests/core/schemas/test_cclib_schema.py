from __future__ import annotations

import gzip
import logging
import os
import shutil
from pathlib import Path
from shutil import copytree

import pytest
from ase.build import bulk
from ase.io import read
from cclib.io import ccread
from maggma.stores import MemoryStore
from monty.json import MontyDecoder, jsanitize

from quacc.calculators.vasp import Vasp
from quacc.schemas.cclib import (
    _cclib_calculate,
    _make_cclib_schema,
    cclib_summarize_run,
)

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True

FILE_DIR = Path(__file__).parent

run1 = FILE_DIR / "test_files" / "gaussian_run1"
log1 = run1 / "Gaussian.log"
run2 = FILE_DIR / "test_files" / "cclib_data"
log2 = run2 / "gau_testopt.log.gz"


@pytest.fixture()
def cclib_obj():
    return ccread(log1)


def setup_module():
    p = FILE_DIR / "test_files" / "cclib_data"

    with (
        gzip.open(p / "psi_test.cube.gz", "r") as f_in,
        open(p / "psi_test.cube", "wb") as f_out,
    ):
        shutil.copyfileobj(f_in, f_out)


def teardown_module():
    p = FILE_DIR / "test_files" / "cclib_data"

    if os.path.exists(p / "psi_test.cube"):
        os.remove(p / "psi_test.cube")


def bad_mock_cclib_calculate(*args, **kwargs):
    msg = "This is a bad run."
    raise ValueError(msg)


def test_cclib_summarize_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    copytree(run1, tmp_path / "test1")
    copytree(run2, tmp_path / "test2")

    # Make sure metadata is made
    atoms = read(log1)
    results = cclib_summarize_run(
        atoms, ".log", directory=tmp_path / "test1", additional_fields={"test": "hi"}
    )
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["spin_multiplicity"] == 1
    assert results["natoms"] == 6
    assert results["attributes"]["metadata"].get("success", None) is True
    assert results["results"].get("energy", None) == pytest.approx(-5516.118738093933)
    assert "pymatgen_version" in results["builder_meta"]

    # Make sure metadata is made
    atoms = read(log2)
    results = cclib_summarize_run(
        atoms, ".log", directory=tmp_path / "test2", additional_fields={"test": "hi"}
    )
    assert results["attributes"]["final_scf_energy"] == pytest.approx(-4091.763)
    assert results["natoms"] == 2
    assert results["charge"] == 0
    assert results["spin_multiplicity"] == 3
    assert results["nelectrons"] == 16
    assert "gau_testopt.log.gz" in results["logfile"]
    assert results.get("attributes") is not None
    assert results["attributes"]["metadata"]["success"] is True
    assert results["input_atoms"]["atoms"][0].position == pytest.approx([0, 0, 0])
    assert results["atoms"][0].position == pytest.approx([0.397382, 0.0, 0.0])
    assert results["attributes"]["homo_energies"] == pytest.approx(
        [-7.054007346511501, -11.618445074798501]
    )
    assert results["attributes"]["lumo_energies"] == pytest.approx(
        [4.2384453353880005, -3.9423854660440005]
    )
    assert results["attributes"]["homo_lumo_gaps"] == pytest.approx(
        [11.292452681899501, 7.6760596087545006]
    )
    assert results["attributes"]["min_homo_lumo_gap"] == pytest.approx(
        7.6760596087545006
    )
    assert len(results["trajectory"]) == 7
    assert results["trajectory"][0] == results["input_atoms"]["atoms"]
    assert results["trajectory"][-1] == results["atoms"]
    assert results["test"] == "hi"

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)

    # Make sure default dir works
    monkeypatch.chdir(tmp_path / "test1")
    cclib_summarize_run(atoms, ".log")

    # Test DB
    atoms = read(log1)
    store = MemoryStore()
    cclib_summarize_run(atoms, ".log", directory=tmp_path / "test1", store=store)
    assert store.count() == 1

    # Make sure info tags are handled appropriately
    atoms = read(log1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = cclib_summarize_run(atoms, ".log", directory=tmp_path / "test1")
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}


def test_errors():
    atoms = bulk("Cu")
    with pytest.raises(ValueError, match="ASE Atoms object has no attached calculator"):
        cclib_summarize_run(atoms, ".log", directory=run1)

    calc = Vasp(atoms)
    atoms.calc = calc
    with pytest.raises(
        ValueError, match="ASE Atoms object's calculator has no results."
    ):
        cclib_summarize_run(atoms, ".log", directory=run1)


def test_cclib_taskdoc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    p = FILE_DIR / "test_files" / "cclib_data"

    # Now we will try two possible extensions, but we will make sure that
    # it fails because the newest log file (.txt) is not valid
    with open(p / "test.txt", "w") as f:
        f.write("I am a dummy log file")
    with pytest.raises(Exception, match="Could not parse") as e:
        doc = _make_cclib_schema(p, [".log", ".txt"])
    os.remove(p / "test.txt")
    assert "Could not parse" in str(e.value)

    # Test a population analysis
    doc = _make_cclib_schema(p, "psi_test.out", analysis="MBO")
    assert doc["pop_analysis"]["mbo"] is not None

    # Let's try with two analysis (also check case-insensitivity)
    doc = _make_cclib_schema(p, "psi_test.out", analysis=["mbo", "density"])
    assert doc["pop_analysis"]["mbo"] is not None
    assert doc["pop_analysis"]["density"] is not None

    # Test a population analysis that will fail
    doc = _make_cclib_schema(p, ".log", analysis="MBO")
    assert doc["pop_analysis"]["mbo"] is None

    doc = _make_cclib_schema(p, "psi_test.out", analysis=["Bader"])
    assert doc["pop_analysis"]["bader"] is not None

    with pytest.raises(FileNotFoundError):
        _make_cclib_schema(p, "does_not_exists.txt")

    # test document can be jsanitized
    d = jsanitize(doc, enum_values=True)

    # and decoded
    MontyDecoder().process_decoded(d)


def test_cclib_calculate(tmp_path, monkeypatch, cclib_obj):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(
        ValueError, match="fakemethod is not a valid cclib population analysis method"
    ):
        _cclib_calculate(cclib_obj, method="fakemethod")

    with pytest.raises(ValueError, match="A cube file must be provided for bader."):
        _cclib_calculate(cclib_obj, method="bader")

    with pytest.raises(FileNotFoundError):
        _cclib_calculate(cclib_obj, method="bader", cube_file="does_not_exists.txt")

    with pytest.raises(FileNotFoundError):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "test_files" / "cclib_data" / "psi_test.cube",
            proatom_dir="does_not_exists",
        )

    with pytest.raises(
        OSError,
        match="PROATOM_DIR environment variable or proatom_dir kwarg needs to be set",
    ):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "test_files" / "cclib_data" / "psi_test.cube",
        )

    with pytest.raises(AssertionError):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "test_files" / "cclib_data" / "psi_test.cube",
            proatom_dir=FILE_DIR / "test_files" / "cclib_data" / "psi_test.cube",
        )


def test_monkeypatches(tmp_path, monkeypatch, cclib_obj, caplog):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(
        "PROATOM_DIR", str(FILE_DIR / "test_files" / "cclib_data" / "proatomdata")
    )
    with pytest.raises(FileNotFoundError):
        _cclib_calculate(
            cclib_obj,
            method="ddec6",
            cube_file=FILE_DIR / "test_files" / "cclib_data" / "psi_test.cube",
        )

    monkeypatch.setattr("cclib.method.Bader.calculate", bad_mock_cclib_calculate)
    with caplog.at_level(logging.WARNING):
        assert (
            _cclib_calculate(
                cclib_obj,
                method="bader",
                cube_file=FILE_DIR / "test_files" / "cclib_data" / "psi_test.cube",
            )
            is None
        )
