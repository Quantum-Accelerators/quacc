from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.mlp.elastic import elastic_tensor_flow

torch = pytest.importorskip("torch")

from importlib.util import find_spec

methods = []
if has_mace := find_spec("mace"):
    methods.append("mace-mp")

if has_matgl := find_spec("matgl"):
    methods.extend(["m3gnet", "chgnet", "tensornet"])

if has_sevennet := find_spec("sevenn"):
    methods.append("sevennet")

if has_orb := find_spec("orb_models"):
    methods.append("orb")

if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if get_token():
        methods.append("fairchem")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.parametrize("method", methods)
def test_elastic_jobs(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_elastic_modulus = {
        "chgnet": 128.184,
        "m3gnet": 126.527,
        "tensornet": 123.915,
        "mace-mp": 130.727,
        "sevennet": 142.296,
        "orb": 190.195,
        "fairchem": 151.367,
    }

    atoms = bulk("Cu")

    outputs = elastic_tensor_flow(
        atoms,
        run_static=False,
        pre_relax=True,
        job_params={
            "all": {"method": method, **calc_kwargs},
            "relax_job": {"opt_params": {"fmax": 0.01}},
        },
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )
    assert outputs["undeformed_result"]["results"]["stress"] == pytest.approx(
        0, abs=1e-2
    )
    assert outputs["elasticity_doc"].bulk_modulus.voigt == pytest.approx(
        ref_elastic_modulus[method], abs=2
    )
    for output in outputs["deformed_results"]:
        assert output["structure_metadata"]["nelements"] == 1
        assert output["structure_metadata"]["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24

    outputs = elastic_tensor_flow(
        atoms,
        run_static=True,
        pre_relax=True,
        job_params={
            "all": {"method": method, **calc_kwargs},
            "relax_job": {"opt_params": {"fmax": 0.01}},
        },
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )

    for output in outputs["deformed_results"]:
        assert output["structure_metadata"]["nelements"] == 1
        assert output["structure_metadata"]["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24

    outputs = elastic_tensor_flow(
        atoms,
        run_static=True,
        pre_relax=False,
        job_params={
            "all": {"method": method, **calc_kwargs},
            "relax_job": {"opt_params": {"fmax": 0.01}},
        },
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )

    for output in outputs["deformed_results"]:
        assert output["structure_metadata"]["nelements"] == 1
        assert output["structure_metadata"]["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24

    outputs = elastic_tensor_flow(
        atoms,
        run_static=False,
        pre_relax=False,
        job_params={
            "all": {"method": method, **calc_kwargs},
            "relax_job": {"opt_params": {"fmax": 0.01}},
        },
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )

    for output in outputs["deformed_results"]:
        assert output["structure_metadata"]["nelements"] == 1
        assert output["structure_metadata"]["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24
