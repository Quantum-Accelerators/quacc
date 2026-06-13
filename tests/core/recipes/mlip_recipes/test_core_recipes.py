from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

import numpy as np
from ase.build import bulk

from quacc.recipes.mlip.core import relax_job, static_job

libraries = []
if has_matcalc := find_spec("matcalc"):
    libraries.append("matcalc")


if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if get_token():
        libraries.append("fairchem")


@pytest.mark.parametrize("library", libraries)
def test_static_job(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    if library == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1p1", "task_name": "omat"}
    elif library == "matcalc":
        calc_kwargs = {"name": "TensorNet-PES-MatPES-PBE-2025.2"}
    else:
        calc_kwargs = {}

    ref_energy = {"matcalc": -3.7303245067596436, "fairchem": -3.7501682869643735}
    atoms = bulk("Cu")
    output = static_job(atoms, library=library, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[library], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.parametrize("library", libraries)
def test_relax_job(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    if library == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1p1", "task_name": "omat"}
    elif library == "matcalc":
        calc_kwargs = {"name": "TensorNet-PES-MatPES-PBE-2025.2"}
    else:
        calc_kwargs = {}

    ref_energy = {"matcalc": -29.842527389526367, "fairchem": -30.001143639922756}

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, library=library, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[library], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.parametrize("library", libraries)
def test_relax_cell_job(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    if library == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1p1", "task_name": "omat"}
    elif library == "matcalc":
        calc_kwargs = {"name": "TensorNet-PES-MatPES-PBE-2025.2"}
    else:
        calc_kwargs = {}

    ref_energy = {"matcalc": -29.87679100036621, "fairchem": -30.005004590392726}

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, library=library, relax_cell=True, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[library], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())


def test_old_imports():
    from quacc.recipes.mlp import _base  # noqa: F401
    from quacc.recipes.mlp.core import relax_job, static_job  # noqa: F401
    from quacc.recipes.mlp.elastic import elastic_tensor_flow  # noqa: F401
    from quacc.recipes.mlp.phonons import phonon_flow  # noqa: F401
