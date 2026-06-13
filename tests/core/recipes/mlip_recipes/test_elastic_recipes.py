from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.mlip.elastic import elastic_tensor_flow

torch = pytest.importorskip("torch")

from importlib.util import find_spec

libraries = []
if has_matcalc := find_spec("matcalc"):
    libraries.append("matcalc")


if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if get_token():
        libraries.append("fairchem")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.parametrize("library", libraries)
def test_elastic_jobs(tmp_path, monkeypatch, library):
    monkeypatch.chdir(tmp_path)

    if library == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1p1", "task_name": "omat"}
    elif library == "matcalc":
        calc_kwargs = {"name": "TensorNet-PES-MatPES-PBE-2025.2"}
    else:
        calc_kwargs = {}

    ref_elastic_modulus = {"matcalc": 132.925, "fairchem": 151.367}

    atoms = bulk("Cu")

    outputs = elastic_tensor_flow(
        atoms,
        run_static=False,
        pre_relax=True,
        job_params={
            "all": {"library": library, **calc_kwargs},
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
        ref_elastic_modulus[library], abs=2
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
            "all": {"library": library, **calc_kwargs},
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
            "all": {"library": library, **calc_kwargs},
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
            "all": {"library": library, **calc_kwargs},
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
