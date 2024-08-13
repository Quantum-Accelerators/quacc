from __future__ import annotations

import pytest

dask = pytest.importorskip("dask")
pytest.importorskip("distributed")

from ase.build import bulk
from dask.distributed import get_client

from quacc import flow, job
from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412

client = get_client()


@pytest.mark.parametrize("job_decorators", [None, {"relax_job": job()}])
def test_functools(tmp_path, monkeypatch, job_decorators):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    delayed = bulk_to_slabs_flow(
        atoms,
        run_static=False,
        job_params={"relax_job": {"opt_params": {"fmax": 0.1}}},
        job_decorators=job_decorators,
    )
    result = client.compute(delayed).result()
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

    assert "atoms" in client.compute(myflow(atoms)).result()


def test_dask_phonon_flow(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    pytest.importorskip("seekpath")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    future = phonon_flow(atoms)
    assert client.compute(future).result()["results"]["thermal_properties"][
        "temperatures"
    ].shape == (101,)


def test_dask_phonon_flow_multistep(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    pytest.importorskip("seekpath")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    future = phonon_flow(relaxed["atoms"])
    assert client.compute(future).result()["results"]["thermal_properties"][
        "temperatures"
    ].shape == (101,)
