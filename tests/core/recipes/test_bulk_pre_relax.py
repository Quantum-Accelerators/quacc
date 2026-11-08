from __future__ import annotations

from functools import wraps
from importlib import import_module

import pytest
from ase.build import bulk

from quacc.recipes.common.elastic import elastic_tensor_flow
from quacc.recipes.common.eos import eos_flow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.wflow_tools.customizers import update_parameters


@pytest.mark.parametrize("workflow", [eos_flow, elastic_tensor_flow])
@pytest.mark.parametrize("run_static", [False, True])
@pytest.mark.parametrize(
    ("flow_params", "job_relax_cell", "cell_relaxed"),
    [
        ({}, False, True),
        ({"job_params": {"relax_job": {"relax_cell": False}}}, True, False),
        ({"pre_relax": False}, True, False),
    ],
)
def test_common_pre_relax(
    tmp_path,
    monkeypatch,
    workflow,
    run_static,
    flow_params,
    job_relax_cell,
    cell_relaxed,
):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu", a=3.5)
    customized_relax = update_parameters(relax_job, {"relax_cell": job_relax_cell})
    job_params = {
        "static_job": {"asap_cutoff": True},
        **flow_params.get("job_params", {}),
    }

    def mark_job(func):
        @wraps(func)
        def wrapped(*args, **kwargs):
            result = func(*args, **kwargs)
            result["custom_decorator"] = True
            return result

        return wrapped

    outputs = workflow(
        atoms,
        customized_relax,
        static_job,
        run_static=run_static,
        **{**flow_params, "job_params": job_params},
        job_decorators={"all": mark_job},
    )

    undeformed = outputs["undeformed_result"]
    assert undeformed["custom_decorator"] is True
    assert (
        undeformed["atoms"].get_volume() != pytest.approx(atoms.get_volume())
    ) is cell_relaxed
    assert undeformed["name"] == (
        "EMT Static"
        if run_static or not flow_params.get("pre_relax", True)
        else "EMT Relax"
    )
    if undeformed["name"] == "EMT Static":
        assert undeformed["parameters"]["asap_cutoff"] is True

    # Pre-relaxation settings must never release the strained cells.
    if workflow is eos_flow:
        from quacc.atoms.deformation import make_volume_strains_from_bulk

        deformed_atoms = make_volume_strains_from_bulk(undeformed["atoms"])
    else:
        deformed_atoms = [
            deformed.to_ase_atoms() for deformed in outputs["deformed_structure_set"]
        ]
    for expected, result in zip(
        deformed_atoms, outputs["deformed_results"]["relax"], strict=True
    ):
        assert result["atoms"].cell.array == pytest.approx(expected.cell.array)
        assert result["custom_decorator"] is True
    if run_static:
        for relaxed, static in zip(
            outputs["deformed_results"]["relax"],
            outputs["deformed_results"]["static"],
            strict=True,
        ):
            assert static["custom_decorator"] is True
            assert static["parameters"]["asap_cutoff"] is True
            assert static["atoms"].cell.array == pytest.approx(
                relaxed["atoms"].cell.array
            )
    else:
        assert outputs["deformed_results"]["static"] == []


@pytest.mark.parametrize("backend", ["emt", "mlip"])
@pytest.mark.parametrize(
    ("job_params", "cell_relaxed"),
    [
        (None, True),
        ({"relax_job": {"relax_cell": False}}, False),
        ({"all": {"relax_cell": True}, "relax_job": {"relax_cell": False}}, False),
        ({"all": {"relax_cell": False}, "relax_job": {"relax_cell": True}}, True),
    ],
)
def test_elastic_wrapper_pre_relax(
    tmp_path, monkeypatch, backend, job_params, cell_relaxed
):
    monkeypatch.chdir(tmp_path)
    recipe = import_module(f"quacc.recipes.{backend}.elastic")
    # Exercise both wrappers without downloading an MLIP model.
    monkeypatch.setattr(recipe, "relax_job", relax_job)
    monkeypatch.setattr(recipe, "static_job", static_job)
    atoms = bulk("Cu", a=3.5)

    outputs = recipe.elastic_tensor_flow(atoms, job_params=job_params)

    assert (
        outputs["undeformed_result"]["atoms"].get_volume()
        != pytest.approx(atoms.get_volume())
    ) is cell_relaxed
    for deformed, result in zip(
        outputs["deformed_structure_set"],
        outputs["deformed_results"]["relax"],
        strict=True,
    ):
        assert result["atoms"].cell.array == pytest.approx(
            deformed.to_ase_atoms().cell.array
        )
