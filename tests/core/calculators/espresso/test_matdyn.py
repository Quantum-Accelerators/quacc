from __future__ import annotations

from copy import deepcopy
from logging import WARNING

import numpy as np
import pytest
from ase import Atoms
from ase.io.espresso import read_fortran_namelist
from numpy.testing import assert_allclose, assert_array_equal

from quacc.calculators.espresso.espresso import Espresso, EspressoTemplate
from quacc.calculators.espresso.utils import read_ifc_has_zstar, read_matdyn_freq

# Header of a QE 7.5 q2r.x force-constant file for NaCl (ibrav=0, polar).
NACL_FC_HEADER = """\
  2    2   0  7.4877033  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
      0.000000000    0.707106781    0.707106781
      0.707106781    0.000000000    0.707106781
      0.707106781    0.707106781    0.000000000
           1  'Na    '    20953.892860652242
           2  'Cl    '    32310.698418201875
    1    1      0.0000000000      0.0000000000      0.0000000000
    2    2      0.7071067811      0.0000000000      0.0000000000
 T   1.4201607893889505
          2.573064743418         -0.000000000000          0.000000000000
"""
# Header of a QE 7.5 force-constant file for Si written with ibrav=2, which
# omits the lattice vectors.
SI_IBRAV2_FC_HEADER = """\
  1    2   2 10.2600000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
           1  'Si    '    25598.823011949724
    1    1      0.0000000000      0.0000000000      0.0000000000
    2    1      0.2500000000      0.2500000000      0.2500000000
 T   2.6664594577970391
         24.799410200025         -0.000000000000          0.000000000000
"""
BAND_POINTS = [
    [0.5, 0.5, 0.5, 20],
    [0, 0, 0, 20],
    [0.5, 0, 0.5, 20],
    [0.375, 0.375, 0.75, 20],
    [0, 0, 0, 1],
]


@pytest.mark.parametrize("native", [False, "list", "string"])
def test_matdyn_band_input(tmp_path, monkeypatch, native):
    """Pass QE endpoint counts through ASE without reinterpreting them."""
    monkeypatch.chdir(tmp_path)
    cards = ["5"] + [" ".join(map(str, point)) for point in BAND_POINTS]
    kwargs = {"qpts": BAND_POINTS}
    if native:
        kwargs = {
            "additional_cards": "\n".join(cards) + "\n" if native == "string" else cards
        }
    original = deepcopy(kwargs)
    calc = Espresso(
        template=EspressoTemplate("matdyn"),
        input_data={
            "input": {
                "asr": "crystal",
                "dos": False,
                "q_in_band_form": True,
                "q_in_cryst_coord": True,
            }
        },
        **kwargs,
    )
    # Repeated input writes must preserve qpts and avoid duplicate cards.
    for _ in range(2):
        calc.write_inputfiles(Atoms(), ["energy"])
        with (tmp_path / "matdyn.in").open() as handle:
            namelist, written = read_fortran_namelist(handle)
        inputs = namelist["input"]
        assert inputs["q_in_band_form"] is True
        assert inputs["q_in_cryst_coord"] is True
        assert inputs["dos"] is False
        assert not {"nk1", "nk2", "nk3"}.intersection(inputs)
        assert written[0] == "5"
        assert written[-1] == "EOF"
        assert_array_equal(np.loadtxt(written[1:-1]), BAND_POINTS)
    assert kwargs == original
    if not native:
        assert calc.parameters["qpts"] == BAND_POINTS


@pytest.mark.parametrize("band_form", [False, True])
def test_matdyn_disconnected_input(tmp_path, monkeypatch, band_form):
    """A zero count skips a segment; explicit lists never add connecting points."""
    monkeypatch.chdir(tmp_path)
    points = [[0, 0, 0, 2], [0.5, 0, 0, 0], [0.5, 0.5, 0, 3], [0, 0, 0, 1]]
    if not band_form:
        points = [[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0, 0]]
    calc = Espresso(
        template=EspressoTemplate("matdyn"),
        input_data={"input": {"q_in_band_form": band_form, "q_in_cryst_coord": True}},
        qpts=np.array(points),
    )
    calc.write_inputfiles(Atoms(), ["energy"])
    with (tmp_path / "matdyn.in").open() as handle:
        _, cards = read_fortran_namelist(handle)
    assert int(cards[0]) == len(points)
    assert_array_equal(np.loadtxt(cards[1:-1]), points)


@pytest.mark.parametrize(
    ("points", "inputs", "extra", "message"),
    [
        ([], {}, {}, "nonempty"),
        ([0, 0, 0], {}, {}, "nonempty"),
        ([[0, 0, 0, 1]], {}, {}, "n, 3"),
        ([[0, 0, 0]], {"q_in_band_form": True}, {}, "n, 4"),
        ([[0, np.nan, 0]], {}, {}, "finite"),
        ([[0, 0, 0, -1]], {"q_in_band_form": True}, {}, "nonnegative integers"),
        ([[0, 0, 0, 1.5]], {"q_in_band_form": True}, {}, "nonnegative integers"),
        ([[0, 0, 0]], {"dos": True}, {}, "dos=True"),
        ([[0, 0, 0]], {"readtau": True}, {}, "additional_cards"),
        ([[0, 0, 0]], {}, {"additional_cards": ["1", "0 0 0"]}, "not both"),
    ],
)
def test_matdyn_invalid_input(tmp_path, monkeypatch, points, inputs, extra, message):
    monkeypatch.chdir(tmp_path)
    calc = Espresso(
        template=EspressoTemplate("matdyn"),
        input_data={"input": inputs},
        qpts=points,
        **extra,
    )
    with pytest.raises(ValueError, match=message):
        calc.write_inputfiles(Atoms(), ["energy"])


@pytest.mark.parametrize("four_columns", [False, True])
@pytest.mark.parametrize("with_gp", [False, True])
def test_matdyn_multiline_frequencies(tmp_path, four_columns, with_gp):
    """Keep point/mode order, negative values, and optional QE distances."""
    weight = " 0.000000" if four_columns else ""
    (tmp_path / "matdyn.freq").write_text(
        " &plot nbnd=9,\n nks=2 /\n"
        f" 0.000000 0.000000 0.000000{weight}\n"
        " -2.0000 3.0000 1.0000 4.0000 5.0000 6.0000\n"
        " 7.0000 8.0000 9.0000\n\n"
        f" 0.500000 0.000000 0.000000{weight}\n"
        " -1.0D+01 9.0d+00 8 7 6 5\n"
        " 4 3 2\n"
    )
    if with_gp:
        (tmp_path / "matdyn.freq.gp").write_text(
            "0 -2 3 1 4 5 6 7 8 9\n0.5 -10 9 8 7 6 5 4 3 2\n"
        )
    results = EspressoTemplate("matdyn").read_results(tmp_path)
    parsed = results["matdyn_results"]
    assert results["energy"] is None
    assert_array_equal(
        parsed["frequencies"],
        [[-2, 3, 1, 4, 5, 6, 7, 8, 9], [-10, 9, 8, 7, 6, 5, 4, 3, 2]],
    )
    assert_array_equal(parsed["qpoints"], [[0, 0, 0], [0.5, 0, 0]])
    assert parsed["frequency_units"] == "cm^-1"
    assert parsed["qpoint_coordinate_system"] == "cartesian"
    assert parsed["qpoint_units"] == "2pi/alat"
    if with_gp:
        assert_array_equal(parsed["distances"], [0, 0.5])
        assert parsed["distance_units"] == "2pi/alat"
    else:
        assert "distances" not in parsed


def test_matdyn_single_point_and_adjacent_negative_values(tmp_path):
    path = tmp_path / "matdyn.freq"
    path.write_text(" &plot nbnd=3, nks=1 /\n 0 0 0\n-1234.0000-2345.0000    3.0000\n")
    parsed = read_matdyn_freq(path)
    assert_array_equal(parsed["frequencies"], [[-1234, -2345, 3]])
    assert parsed["qpoints"].shape == (1, 3)


@pytest.mark.parametrize("custom_names", [False, True])
def test_matdyn_dos_and_filenames(tmp_path, monkeypatch, custom_names):
    monkeypatch.chdir(tmp_path)
    freq_name = "Custom.freq" if custom_names else "matdyn.freq"
    dos_name = "Custom.dos" if custom_names else "matdyn.dos"
    inputs = {"dos": True, "nk1": 4, "nk2": 4, "nk3": 4}
    if custom_names:
        inputs.update(flfrq=freq_name, fldos=dos_name)
        # Ignore unrelated files left in the directory.
        (tmp_path / "matdyn.freq").write_text("not the requested frequencies")
        (tmp_path / "matdyn.dos").write_text("not the requested DOS")
    calc = Espresso(template=EspressoTemplate("matdyn"), input_data={"input": inputs})
    calc.write_inputfiles(Atoms(), ["energy"])
    with (tmp_path / "matdyn.in").open() as handle:
        written, cards = read_fortran_namelist(handle)
    assert written["input"]["flfrq"] == freq_name
    assert written["input"]["fldos"] == dos_name
    assert all(written["input"][key] == 4 for key in ("nk1", "nk2", "nk3"))
    assert cards == ["EOF"]
    (tmp_path / dos_name).write_text("# frequency dos integrated\n0 1 2\n1 2 3\n")
    (tmp_path / "q2r.fc").write_text(NACL_FC_HEADER)
    # Reading a DOS-only calculation keeps the original results layout.
    parsed = EspressoTemplate("matdyn").read_results(tmp_path)["matdyn_results"]
    assert set(parsed) == {"phonon_dos"}
    assert_array_equal(parsed["phonon_dos"], [[0, 1, 2], [1, 2, 3]])
    (tmp_path / freq_name).write_text("DOS frequency files need not be parsed")
    parsed = EspressoTemplate("matdyn").read_results(tmp_path)["matdyn_results"]
    assert set(parsed) == {"phonon_dos"}


@pytest.mark.parametrize("freq_name", ["matdyn.freq", "Custom.freq"])
def test_matdyn_custom_dispersion_filename(tmp_path, monkeypatch, freq_name):
    monkeypatch.chdir(tmp_path)
    calc = Espresso(
        template=EspressoTemplate("matdyn"),
        input_data={"input": {"dos": False, "flfrq": freq_name}},
        qpts=[[0, 0, 0]],
    )
    calc.write_inputfiles(Atoms(), ["energy"])
    if freq_name != "matdyn.freq":
        (tmp_path / "matdyn.freq").write_text("ignore default when using a custom name")
    (tmp_path / freq_name).write_text(" &plot nbnd=3, nks=1 /\n0 0 0\n-1 2 3\n")
    (tmp_path / (freq_name + ".gp")).write_text("0 -1 2 3\n")
    parsed = EspressoTemplate("matdyn").read_results(tmp_path)["matdyn_results"]
    assert_array_equal(parsed["frequencies"], [[-1, 2, 3]])
    assert_array_equal(parsed["distances"], [0])


def test_matdyn_gp_only(tmp_path):
    (tmp_path / "matdyn.freq.gp").write_text("# distance frequencies\n0 -1 2 3\n")
    parsed = read_matdyn_freq(tmp_path / "matdyn.freq")
    assert_array_equal(parsed["frequencies"], [[-1, 2, 3]])
    assert_array_equal(parsed["distances"], [0])
    assert "qpoints" not in parsed


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("0 0 0\n1 2 3\n", "header"),
        (" &plot nbnd=0, nks=1 /\n", "dimensions"),
        (" &plot nbnd=3, nks=2 /\n0 0 0\n1 2 3\n", "Incomplete"),
        (" &plot nbnd=3, nks=1 /\n0 0\n1 2 3\n", "coordinates"),
        (" &plot nbnd=3, nks=1 /\n0 0 0\n1 2 3 4\n", "number"),
        (" &plot nbnd=3, nks=1 /\n0 0 0\n1 2 3\n0 0 0\n", "extra"),
        (" &plot nbnd=3, nks=1 /\n0 0 0\n********** 2 3\n", "numeric"),
    ],
)
def test_matdyn_malformed_freq(tmp_path, content, message):
    path = tmp_path / "matdyn.freq"
    path.write_text(content)
    with pytest.raises(ValueError, match=message):
        read_matdyn_freq(path)


@pytest.mark.parametrize("gp_content", ["0 1 2\n", "0 3 2 1\n", "0 1 2 3\n1 1 2 3\n"])
def test_matdyn_inconsistent_gp(tmp_path, gp_content):
    path = tmp_path / "matdyn.freq"
    path.write_text(" &plot nbnd=3, nks=1 /\n0 0 0\n1 2 3\n")
    path.with_name("matdyn.freq.gp").write_text(gp_content)
    with pytest.raises(ValueError, match="disagree"):
        read_matdyn_freq(path)


def test_matdyn_gp_disconnected_order(tmp_path):
    """Distances in QE's .gp include jumps; do not infer segment connectivity."""
    (tmp_path / "matdyn.freq").write_text(
        " &plot nbnd=3, nks=4 /\n"
        "0 0 0 0\n-1 2 3\n0.5 0 0 0.5\n4 5 6\n"
        "0 0.5 0 0.5\n7 8 9\n0 0.75 0 0.75\n10 11 12\n"
    )
    (tmp_path / "matdyn.freq.gp").write_text(
        "0 -1 2 3\n0.5 4 5 6\n1.207107 7 8 9\n1.457107 10 11 12\n"
    )
    parsed = read_matdyn_freq(tmp_path / "matdyn.freq")
    assert_array_equal(
        parsed["qpoints"], [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0.75, 0]]
    )
    assert_allclose(parsed["distances"], [0, 0.5, 1.207107, 1.457107])
    assert_array_equal(parsed["frequencies"][:, 0], [-1, 4, 7, 10])


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (NACL_FC_HEADER, True),
        (NACL_FC_HEADER.replace(" T   1.42", " F   1.42"), False),
        (NACL_FC_HEADER.replace(" T   1.42", " .TRUE.   1.42"), True),
        (SI_IBRAV2_FC_HEADER, True),
        (SI_IBRAV2_FC_HEADER.replace(" T   2.66", " F   2.66"), False),
    ],
)
def test_read_ifc_has_zstar(tmp_path, header, expected):
    (tmp_path / "q2r.fc").write_text(header)
    assert read_ifc_has_zstar(tmp_path / "q2r.fc") is expected


def test_read_ifc_has_zstar_unknown(tmp_path):
    """XML and missing force-constant files leave the flag unknown."""
    (tmp_path / "q2r.fc.xml").write_text(NACL_FC_HEADER)
    assert read_ifc_has_zstar(tmp_path / "q2r.fc.xml") is None
    assert read_ifc_has_zstar(tmp_path / "missing.fc") is None


@pytest.mark.parametrize(
    "header",
    [
        "",
        "force constants\n",
        NACL_FC_HEADER.replace(" T   1.42", " 1.42"),
        "\n".join(NACL_FC_HEADER.splitlines()[:8]),
    ],
)
def test_read_ifc_has_zstar_unrecognized(tmp_path, caplog, header):
    """matdyn.x already read the file, so an unknown flag must not fail the job."""
    (tmp_path / "q2r.fc").write_text(header)
    with caplog.at_level(WARNING):
        assert read_ifc_has_zstar(tmp_path / "q2r.fc") is None
    assert "Could not read has_zstar" in caplog.text


@pytest.mark.parametrize("with_fc", [False, True])
def test_matdyn_dispersion_has_zstar(tmp_path, monkeypatch, with_fc):
    """Dispersion results report the dielectric flag of the force constants."""
    monkeypatch.chdir(tmp_path)
    calc = Espresso(
        template=EspressoTemplate("matdyn"),
        input_data={"input": {"dos": False}},
        qpts=[[0, 0, 0]],
    )
    calc.write_inputfiles(Atoms(), ["energy"])
    (tmp_path / "matdyn.freq").write_text(" &plot nbnd=3, nks=1 /\n0 0 0\n0 0 0\n")
    if with_fc:
        (tmp_path / "q2r.fc").write_text(NACL_FC_HEADER)
    parsed = EspressoTemplate("matdyn").read_results(tmp_path)["matdyn_results"]
    if with_fc:
        assert parsed["has_zstar"] is True
    else:
        assert "has_zstar" not in parsed
