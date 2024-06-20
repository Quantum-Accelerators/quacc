from __future__ import annotations

import pytest

parsl = pytest.importorskip("parsl")

import os

from ase.build import bulk

from quacc import flow, job
from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


@pytest.mark.parametrize("job_decorators", [None, {"relax_job": job()}])
def test_functools(tmp_path, monkeypatch, job_decorators):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    result = bulk_to_slabs_flow(
        atoms,
        run_static=False,
        job_params={"relax_job": {"opt_params": {"fmax": 0.1}}},
        job_decorators=job_decorators,
    ).result()
    assert len(result) == 4
    assert "atoms" in result[-1]
    assert result[-1]["parameters_opt"]["fmax"] == 0.1


def test_copy_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    @flow
    def myflow(atoms):
        result1 = relax_job(atoms)
        return relax_job(result1["atoms"], copy_files={result1["dir_name"]: "opt.*"})

    assert "atoms" in myflow(atoms).result()


def test_settings_swap(tmp_path_factory):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")
    atoms = bulk("Cu")

    future1 = bulk_to_slabs_flow(
        atoms,
        job_decorators={"relax_job": job(settings_swap={"RESULTS_DIR": tmp_dir1})},
    )
    future2 = bulk_to_slabs_flow(
        atoms,
        job_decorators={"relax_job": job(settings_swap={"RESULTS_DIR": tmp_dir1})},
    )
    future3 = bulk_to_slabs_flow(
        atoms,
        job_decorators={"relax_job": job(settings_swap={"RESULTS_DIR": tmp_dir1})},
    )
    future1.result(), future2.result(), future3.result()

    assert len(os.listdir(tmp_dir1)) == 12


def test_settings_swap_all(tmp_path_factory):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")

    atoms = bulk("Cu")

    future1 = bulk_to_slabs_flow(
        atoms, job_decorators={"all": job(settings_swap={"RESULTS_DIR": tmp_dir1})}
    )
    future2 = bulk_to_slabs_flow(
        atoms, job_decorators={"all": job(settings_swap={"RESULTS_DIR": tmp_dir1})}
    )
    future3 = bulk_to_slabs_flow(
        atoms, job_decorators={"all": job(settings_swap={"RESULTS_DIR": tmp_dir1})}
    )
    future1.result(), future2.result(), future3.result()

    assert len(os.listdir(tmp_dir1)) == 24


def test_phonon_flow(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms)
    assert output.result()["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )


def test_phonon_flow_multistep(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    output = phonon_flow(relaxed["atoms"])
    assert output.result()["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )
