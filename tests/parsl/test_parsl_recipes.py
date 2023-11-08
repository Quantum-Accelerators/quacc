import contextlib

import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job


def setup_module():
    with contextlib.suppress(Exception):
        parsl.load()


def teardown_module():
    parsl.clear()


def test_phonon_flow(tmpdir):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    tmpdir.chdir()
    atoms = bulk("Cu")
    output = phonon_flow(atoms)
    assert output.result()["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )


def test_phonon_flow_multistep(tmpdir):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    tmpdir.chdir()
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    output = phonon_flow(relaxed["atoms"])
    assert output.result()["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )
