import pytest

parsl = pytest.importorskip("parsl")

from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412

# from quacc import flow


DEFAULT_SETTINGS = SETTINGS.model_copy()


@pytest.mark.parametrize("chdir", [True, False])
def test_functools(tmp_path, monkeypatch, chdir):
    monkeypatch.chdir(tmp_path)

    SETTINGS.CHDIR = chdir

    atoms = bulk("Cu")
    result = bulk_to_slabs_flow(
        atoms, job_params={"relax_job": {"opt_params": {"fmax": 0.1}}}, run_static=False
    ).result()
    assert len(result) == 4
    assert "atoms" in result[-1]
    assert result[-1]["fmax"] == 0.1

    SETTINGS.CHDIR = DEFAULT_SETTINGS.CHDIR


# def test_copy_files(tmp_path, monkeypatch):
#     monkeypatch.chdir(tmp_path)
#     atoms = bulk("Cu")

#     @flow
#     def myflow(atoms):
#         result1 = relax_job(atoms)
#         return relax_job(result1["atoms"], copy_files={result1["dir_name"]: "opt.*"})

#     assert "atoms" in myflow(atoms).result()


@pytest.mark.parametrize("chdir", [True, False])
def test_phonon_flow(tmp_path, monkeypatch, chdir):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    SETTINGS.CHDIR = chdir

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms)
    assert output.result()["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )

    SETTINGS.CHDIR = DEFAULT_SETTINGS.CHDIR


@pytest.mark.parametrize("chdir", [True, False])
def test_phonon_flow_multistep(tmp_path, monkeypatch, chdir):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    SETTINGS.CHDIR = chdir

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    output = phonon_flow(relaxed["atoms"])
    assert output.result()["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )

    SETTINGS.CHDIR = DEFAULT_SETTINGS.CHDIR
