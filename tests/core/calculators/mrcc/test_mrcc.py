import re

import numpy as np
import pytest

import shutil

from ase.atoms import Atoms
from quacc.calculators.mrcc.mrcc import get_version_from_mrcc_header, MrccProfile, MRCC

def test_mrcc_version_from_string():
    reference_outputfile = """                              www.mrcc.hu
 
                     Release date: August 28, 2023

 ************************ 2024-05-28 17:14:24 *************************
 """
    version = get_version_from_mrcc_header(reference_outputfile)
    assert version == 'August 28, 2023'


def test_mrcc_version_from_executable():
    # only check the format to be compatible with future versions
    version_regexp = re.compile(r"(January|February|March|April|May|June|July|August|September|October|November|December)\s\d{1,2},\s\d{4}")

    dmrcc_path = shutil.which("dmrcc") 

    MyMrccProfile = MrccProfile(dmrcc_path)

    version = MyMrccProfile.version()

    assert version_regexp.match(version)


def test_mrcc_singlepoint(tmpdir):
    dmrcc_path = shutil.which("dmrcc") 

    MyMrccProfile = MrccProfile(dmrcc_path)

    calc = MRCC(label='mrcc',
                profile=MyMrccProfile,
                mrccinput={'calc':'PBE','basis':'STO-3G'}, mrccblocks="""symm=off""",
                directory=tmpdir)

    #Geometry input. Either like this:
    water = Atoms('H2O',positions=[[1.0,0.0,0.0],[2.0,0.0,0.0],[3.0,0.0,0.0]],calculator=calc)
    energy = water.get_potential_energy()

    assert energy == pytest.approx(-2026.1497783941234)

