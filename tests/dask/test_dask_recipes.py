import gzip
import os
from datetime import datetime
from pathlib import Path

import pytest
from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

dask = pytest.importorskip("dask")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "dask",
    reason="This test requires the Dask workflow engine",
)

from dask.distributed import default_client

client = default_client()


def test_dask_speed(tmp_path, monkeypatch):
    """This test is critical for making sure we are using multiple cores"""
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    delayed = bulk_to_slabs_flow(
        atoms,
        slab_relax_kwargs={
            "opt_params": {"optimizer_kwargs": {"logfile": "test_dask_speed.log"}}
        },
    )
    result = client.gather(client.compute(delayed))
    assert len(result) == 4
    assert "atoms" in result[-1]

    times = []
    fs = os.listdir(SETTINGS.RESULTS_DIR)
    fs.sort()
    assert fs

    for d in fs:
        p = Path(SETTINGS.RESULTS_DIR / d, "test_dask_speed.log.gz")
        if p.exists():
            with gzip.open(p, "rt") as file:
                time = []
                for line in file:
                    if ":" in line:
                        time_format = "%H:%M:%S"
                        time_object = datetime.strptime(line.split()[2], time_format)
                        time.append(time_object)
            times.append(time)

    assert times[1][0] < times[0][-1]


def test_dask_phonon_flow(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    future = phonon_flow(atoms)
    assert client.compute(future).result()["results"]["thermal_properties"][
        "temperatures"
    ].shape == (101,)


def test_dask_phonon_flow_multistep(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    future = phonon_flow(relaxed["atoms"])
    assert client.compute(future).result()["results"]["thermal_properties"][
        "temperatures"
    ].shape == (101,)
