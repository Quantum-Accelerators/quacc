import gzip
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from ase.io.espresso import construct_namelist, read_fortran_namelist

from quacc.calculators.espresso.io import (
    read_espresso_ph,
    write_espresso_io,
    write_espresso_ph,
)
from quacc.calculators.espresso.keys import ALL_KEYS

RUN_PATH = Path(__file__).parent / "test_runs"


def _read_gzip_file(file_path):
    with gzip.open(file_path, "rt", encoding="utf-8") as file:
        content = file.read()
    return content


def test_write_espresso_io():
    input_data = {
        "degauss": 2,
        "Emin": True,
        "Emax": None,
        "nbnd": 10,
        "ecutwfc": 20,
        "irmin(2, 6)": 0.1,
        "outdir(5)": "/path/to/outdir",
    }

    input_data = construct_namelist(input_data, ALL_KEYS["projwfc"])

    string_io = StringIO()

    write_espresso_io(string_io, input_data=input_data)

    expected = (
        "&PROJWFC\n"
        "   outdir(5)        = '/path/to/outdir'\n"
        "   degauss          = 2\n"
        "   emin             = .true.\n"
        "   irmin(2, 6)      = 0.1\n"
        "/\n"
        "EOF"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


def test_write_espresso_ph_single():
    input_data = {
        "amass(1)": 1.0,
        "amass(2)": 2.0,
        "prefix": "prefix",
        "outdir": "/path/to/outdir",
        "eth_rps": 0.1,
    }

    qpts = (0.5, -0.1, 1 / 3)

    input_data = construct_namelist(input_data, ALL_KEYS["ph"])

    string_io = StringIO()

    write_espresso_ph(string_io, input_data=input_data, qpts=qpts)

    expected = (
        "&INPUTPH\n"
        "   amass(1)         = 1.0\n"
        "   amass(2)         = 2.0\n"
        "   outdir           = '/path/to/outdir'\n"
        "   prefix           = 'prefix'\n"
        "   eth_rps          = 0.1\n"
        "/\n"
        "0.50000000 -0.10000000 0.33333333\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


def test_write_espresso_ph_list():
    input_data = {
        "amass(1)": 1.0,
        "amass(2)": 2.0,
        "prefix": "prefix",
        "outdir": "/path/to/outdir",
        "eth_rps": 0.1,
        "qplot": True,
        "ldisp": True,
    }

    qpts = [(0.5, -0.1, 1 / 3, 2), (0.1, 0.2, 0.3, 10), (0.2, 0.3, 0.4, 1)]

    input_data = construct_namelist(input_data, ALL_KEYS["ph"])

    string_io = StringIO()

    write_espresso_ph(string_io, input_data=input_data, qpts=qpts)

    expected = (
        "&INPUTPH\n"
        "   amass(1)         = 1.0\n"
        "   amass(2)         = 2.0\n"
        "   outdir           = '/path/to/outdir'\n"
        "   prefix           = 'prefix'\n"
        "   eth_rps          = 0.1\n"
        "   qplot            = .true.\n"
        "   ldisp            = .true.\n"
        "/\n"
        "3\n"
        "0.50000000 -0.10000000 0.33333333 2\n"
        "0.10000000 0.20000000 0.30000000 10\n"
        "0.20000000 0.30000000 0.40000000 1\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


def test_write_espresso_ph_nat_todo():
    input_data = {
        "amass(1)": 1.0,
        "amass(2)": 2.0,
        "prefix": "prefix",
        "outdir": "/path/to/outdir",
        "eth_rps": 0.1,
        "qplot": True,
        "nat_todo": True,
        "ldisp": True,
    }

    qpts = [(0.5, -0.1, 1 / 3, 1), (0.1, 0.2, 0.3, -1), (0.2, 0.3, 0.4, 4)]

    input_data = construct_namelist(input_data, ALL_KEYS["ph"])

    string_io = StringIO()

    write_espresso_ph(string_io, input_data=input_data, qpts=qpts, nat_todo=[1, 2, 3])

    expected = (
        "&INPUTPH\n"
        "   amass(1)         = 1.0\n"
        "   amass(2)         = 2.0\n"
        "   outdir           = '/path/to/outdir'\n"
        "   prefix           = 'prefix'\n"
        "   eth_rps          = 0.1\n"
        "   qplot            = .true.\n"
        "   ldisp            = .true.\n"
        "   nat_todo         = .true.\n"
        "/\n"
        "3\n"
        "0.50000000 -0.10000000 0.33333333 1\n"
        "0.10000000 0.20000000 0.30000000 -1\n"
        "0.20000000 0.30000000 0.40000000 4\n"
        "1 2 3\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


@pytest.mark.parametrize("i", range(1, 66, 5))
def test_read_espresso_ph_all(i):
    file_contents = _read_gzip_file(RUN_PATH / f"run{i}.txt.gz")
    fd = StringIO(file_contents)
    read_espresso_ph(fd)


def test_read_espresso_ph_1():
    file_contents = _read_gzip_file(RUN_PATH / "run1.txt.gz")
    fd = StringIO(file_contents)
    results = read_espresso_ph(fd)

    assert len(results) == 8
    assert (0, 0, 0) in results
    assert np.unique(results[(0, 0, 0)]["freqs"]).shape[0] == 1
    assert np.unique(results[(0, 0, 0)]["freqs"])[0] == 0.173268
    assert len(results[(0, 0, 0)]["eqpoints"]) == 1
    assert results[(0, 0, 0)]["atoms"].symbols == ["Al"]

    assert (0.75, -0.25, 0.75) in results
    assert np.unique(results[(0.75, -0.25, 0.75)]["freqs"]).shape[0] == 3
    assert np.unique(results[(0.75, -0.25, 0.75)]["freqs"])[2] == 8.791383
    assert len(results[(0.75, -0.25, 0.75)]["eqpoints"]) == 24
    assert results[(0.75, -0.25, 0.75)]["atoms"].symbols == ["Al"]
