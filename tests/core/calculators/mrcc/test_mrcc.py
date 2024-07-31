from __future__ import annotations

import pytest
from ase.atoms import Atoms

from quacc import get_settings
from quacc.calculators.mrcc.mrcc import MRCC, MrccProfile, _get_version_from_mrcc_header


def test_mrcc_version_from_string():
    reference_outputfile = """                              www.mrcc.hu

                     Release date: August 28, 2023

 ************************ 2024-05-28 17:14:24 *************************
 """
    version = _get_version_from_mrcc_header(reference_outputfile)
    assert version == "August 28, 2023"


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
