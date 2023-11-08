import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job

pytest.importorskip("phonopy")

from quacc.recipes.emt.phonons import phonon_flow


def test_phonon_flow(tmpdir):
    tmpdir.chdir()
    atoms = bulk("Cu")
    output = phonon_flow(atoms)
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)


def test_phonon_flow_multistep(tmpdir):
    tmpdir.chdir()
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    output = phonon_flow(relaxed["atoms"])
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)
