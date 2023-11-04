import pytest
from ase.build import bulk

parsl = pytest.importorskip("parsl")


def setup_module():
    import contextlib

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
