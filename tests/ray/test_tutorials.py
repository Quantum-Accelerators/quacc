from __future__ import annotations

import pytest

ray = pytest.importorskip("ray")

from ase.build import bulk, molecule

from quacc.recipes.emt.core import relax_job, static_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


def test_tutorial1a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    future = relax_job(atoms)

    assert "atoms" in future.result()


def test_tutorial1b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    future = bulk_to_slabs_flow(atoms)

    assert "atoms" in future.result()[0]


def test_tutorial2a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def workflow(atoms):
        future1 = relax_job(atoms)
        return static_job(future1.result()["atoms"])

    atoms = bulk("Cu")

    future = workflow(atoms)

    assert "atoms" in future.result()


def test_tutorial2b(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def workflow(atoms1, atoms2):
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)
        return [result1, result2]

    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    futures = workflow(atoms1, atoms2)

    results = [f.result() for f in futures]

    assert "atoms" in results[0]
    assert "atoms" in results[1]


def test_tutorial2c(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def workflow(atoms):
        relaxed_bulk = (relax_job(atoms)).result()
        return bulk_to_slabs_flow(relaxed_bulk["atoms"], run_static=False)

    atoms = bulk("Cu")

    future = workflow(atoms)

    result = future.result()

    assert len(result) == 4
    assert "atoms" in result[0]
