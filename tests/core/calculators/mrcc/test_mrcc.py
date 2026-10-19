from __future__ import annotations

import pytest
from ase.atoms import Atoms

from quacc import get_settings
from quacc.calculators.mrcc import mrcc
from quacc.calculators.mrcc.mrcc import (
    MRCC,
    MrccProfile,
    MrccTemplate,
    _get_version_from_mrcc_header,
)

MRCC_EXECUTE = MrccTemplate.execute


def test_mrcc_version_from_string():
    reference_outputfile = """                              www.mrcc.hu

                     Release date: August 28, 2023

 ************************ 2024-05-28 17:14:24 *************************
 """
    version = _get_version_from_mrcc_header(reference_outputfile)
    assert version == "August 28, 2023"


def test_mrcc_profile(monkeypatch):
    header = "Release date: August 28, 2023\n"
    monkeypatch.setattr(mrcc, "read_stdout", lambda command: header)
    profile = MrccProfile(command="dmrcc")

    assert profile.version() == "August 28, 2023"
    assert profile.get_calculator_command("MINP") == ["MINP"]


def test_mrcc_template_execute(tmp_path):
    class Profile:
        def run(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    profile = Profile()
    MRCC_EXECUTE(MrccTemplate(), tmp_path, profile)

    assert profile.args == (tmp_path, "MINP", "mrcc.out")
    assert profile.kwargs == {"errorfile": "mrcc.err"}


def test_mrcc_template_load_profile(monkeypatch):
    expected = object()
    monkeypatch.setattr(
        MrccProfile, "from_config", lambda cfg, name, **kwargs: expected
    )

    assert MrccTemplate().load_profile(object()) is expected


def test_mrcc_singlepoint(tmp_path):
    calc = MRCC(
        profile=MrccProfile(command=get_settings().MRCC_CMD),
        calc="PBE",
        basis="STO-3G",
        symm="off",
        directory=tmp_path,
    )

    # Geometry input. Either like this:
    water = Atoms(
        "H2O",
        positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        calculator=calc,
    )
    energy = water.get_potential_energy()

    assert energy == pytest.approx(-2026.1497783941234)
