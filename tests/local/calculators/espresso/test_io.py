from __future__ import annotations

from io import StringIO
from pathlib import Path

from ase.io.espresso import (
    construct_namelist,
    read_fortran_namelist,
)

from quacc.calculators.espresso.io import write_espresso_io, write_espresso_ph
from quacc.calculators.espresso.keys import ALL_KEYS


def test_write_espresso_io():
    input_data = {
        "degauss": 2,
        "Emin": True,
        "Emax": None,
        "nbnd": 10,
        "ecutwfc": 20,
        "irmin(2, 6)": 0.1,
        "outdir(5)": str(Path("/path/to/outdir")),
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
        "outdir": str(Path("/path/to/outdir")),
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
        "outdir": str(Path("/path/to/outdir")),
        "eth_rps": 0.1,
        "qplot": True,
    }

    qpts = [(0.5, -0.1, 1 / 3), (0.1, 0.2, 0.3), (0.2, 0.3, 0.4)]

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
        "/\n"
        "3\n"
        "0.50000000 -0.10000000 0.33333333\n"
        "0.10000000 0.20000000 0.30000000\n"
        "0.20000000 0.30000000 0.40000000\n"
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
        "outdir": str(Path("/path/to/outdir")),
        "eth_rps": 0.1,
        "qplot": True,
        "nat_todo": True,
    }

    qpts = [(0.5, -0.1, 1 / 3), (0.1, 0.2, 0.3), (0.2, 0.3, 0.4)]

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
        "   nat_todo         = .true.\n"
        "/\n"
        "3\n"
        "0.50000000 -0.10000000 0.33333333\n"
        "0.10000000 0.20000000 0.30000000\n"
        "0.20000000 0.30000000 0.40000000\n"
        "1 2 3\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected
