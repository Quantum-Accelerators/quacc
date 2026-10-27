from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk

from quacc import flow
from quacc.recipes.common.eos import eos_flow, eos_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.wflow_tools.customizers import update_parameters


def _get_fit(eos_doc, model):
    return next(fit for fit in eos_doc.eos if fit.model.value == model)


def test_eos_subflow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    relax_result = relax_job(atoms, relax_cell=True)
    static_result = static_job(relax_result["atoms"])
    v_static = static_result["atoms"].get_volume()
    volume_strains = [-0.075, -0.05, -0.025, 0.025, 0.05, 0.075]

    outputs = eos_subflow(
        static_result,
        static_job=static_job,
        strain_kwargs={"volume_strains": volume_strains},
    )

    # The input result is passed through unchanged
    assert outputs["undeformed_result"] is static_result

    # One static calculation per strain, on the strained volume
    assert outputs["deformed_results"]["relax"] == []
    statics = outputs["deformed_results"]["static"]
    assert len(statics) == len(volume_strains)
    for strain, output in zip(volume_strains, statics, strict=True):
        assert output["name"] == "EMT Static"
        assert output["atoms"].get_volume() == pytest.approx((1 + strain) * v_static)

    # With include_undeformed=True (default), the static result is an extra fit point
    fit_results = [static_result, *statics]
    expected_volumes = sorted(r["atoms"].get_volume() for r in fit_results)
    assert outputs["eos_doc"].volumes == pytest.approx(expected_volumes)

    # Regression value, close to the elastic recipe's 134.579 GPa for EMT Cu, and
    # the fitted minimum sits at the relaxed volume
    fit = _get_fit(outputs["eos_doc"], "birch_murnaghan")
    assert pytest.approx(134.299, abs=1e-2) == fit.B0
    assert pytest.approx(v_static, rel=1e-3) == fit.V0


def test_eos_subflow_from_relax(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    relax_result = relax_job(atoms, relax_cell=True)
    v_relax = relax_result["atoms"].get_volume()
    volume_strains = [-0.05, -0.025, 0.025, 0.05]

    outputs = eos_subflow(
        relax_result,
        static_job=update_parameters(static_job, {"asap_cutoff": True}),
        strain_kwargs={"volume_strains": volume_strains},
        include_undeformed=False,
    )

    assert outputs["undeformed_result"] is relax_result
    statics = outputs["deformed_results"]["static"]
    for strain, output in zip(volume_strains, statics, strict=True):
        # The customized static_job is used for every strained structure
        assert output["parameters"]["asap_cutoff"] is True
        assert output["atoms"].get_volume() == pytest.approx((1 + strain) * v_relax)

    # With include_undeformed=False, only the strained statics are fit
    expected_volumes = [(1 + strain) * v_relax for strain in volume_strains]
    assert outputs["eos_doc"].volumes == pytest.approx(expected_volumes)


def test_eos_subflow_relax(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # A rattled cell, so the fixed-cell relaxations have atoms to move
    atoms = bulk("Cu", cubic=True)
    atoms.rattle(stdev=0.05, seed=42)
    static_result = static_job(atoms)
    v0 = atoms.get_volume()
    volume_strains = [-0.05, -0.025, 0.025, 0.05]

    outputs = eos_subflow(
        static_result,
        relax_job=relax_job,
        static_job=static_job,
        strain_kwargs={"volume_strains": volume_strains},
    )

    relaxes = outputs["deformed_results"]["relax"]
    statics = outputs["deformed_results"]["static"]
    assert len(relaxes) == len(statics) == len(volume_strains)
    for strain, relax, static in zip(volume_strains, relaxes, statics, strict=True):
        assert relax["name"] == "EMT Relax"
        assert static["name"] == "EMT Static"

        # Only the atomic positions relax; the strained cell is kept
        strained_cell = atoms.cell.array * (1 + strain) ** (1 / 3)
        assert relax["atoms"].cell.array == pytest.approx(strained_cell)
        assert relax["atoms"].get_volume() == pytest.approx((1 + strain) * v0)
        assert not np.allclose(
            relax["atoms"].get_scaled_positions(), atoms.get_scaled_positions()
        )

        # The static runs on the relaxed structure
        assert static["atoms"].positions == pytest.approx(relax["atoms"].positions)

    # The fit uses the static energies (plus the undeformed result)
    fit_energies = dict(
        zip(outputs["eos_doc"].volumes, outputs["eos_doc"].energies, strict=True)
    )
    for static in statics:
        volume = static["atoms"].get_volume()
        energy = next(e for v, e in fit_energies.items() if v == pytest.approx(volume))
        assert energy == pytest.approx(static["results"]["energy"])


def test_eos_subflow_relax_only(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    static_result = static_job(bulk("Cu"))
    volume_strains = [-0.05, -0.025, 0.025, 0.05]

    outputs = eos_subflow(
        static_result,
        relax_job=relax_job,
        strain_kwargs={"volume_strains": volume_strains},
        include_undeformed=False,
    )

    relaxes = outputs["deformed_results"]["relax"]
    assert outputs["deformed_results"]["static"] == []
    assert len(relaxes) == len(volume_strains)

    # Without a static job, the relaxation energies are fit
    assert outputs["eos_doc"].energies == pytest.approx(
        [r["results"]["energy"] for r in relaxes]
    )


def test_eos_subflow_bad_inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    static_result = static_job(bulk("Cu"))

    # 3 strains plus the undeformed result is the minimum of 4 points
    outputs = eos_subflow(
        static_result,
        static_job=static_job,
        strain_kwargs={"volume_strains": [-0.05, 0.025, 0.05]},
    )
    assert len(outputs["eos_doc"].volumes) == 4

    with pytest.raises(ValueError, match="At least 4"):
        eos_subflow(
            static_result,
            static_job=static_job,
            strain_kwargs={"volume_strains": [-0.05, 0.025, 0.05]},
            include_undeformed=False,
        )

    with pytest.raises(ValueError, match="relax_job or static_job"):
        eos_subflow(static_result)


def test_eos_subflow_in_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # The subflow receives the static result from earlier jobs in the same flow
    @flow
    def workflow(atoms):
        relax_result = relax_job(atoms, relax_cell=True)
        static_result = static_job(relax_result["atoms"])
        return eos_subflow(static_result, static_job=static_job)

    outputs = workflow(bulk("Cu"))
    assert len(outputs["deformed_results"]["static"]) == 10
    assert len(outputs["eos_doc"].volumes) == 11
    fit = _get_fit(outputs["eos_doc"], "birch_murnaghan")
    assert pytest.approx(134.311, abs=1e-2) == fit.B0


def test_eos_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Start from a compressed cell so the pre-relaxation has to change the volume
    atoms = bulk("Cu", a=3.5)
    outputs = eos_flow(atoms, relax_job, static_job)

    # Default: a cell relaxation on the bulk structure
    undeformed = outputs["undeformed_result"]
    assert undeformed["name"] == "EMT Relax"
    assert undeformed["atoms"].get_volume() != pytest.approx(atoms.get_volume())

    # Default: only fixed-cell relaxations on the strained structures
    assert len(outputs["deformed_results"]["relax"]) == 10
    assert outputs["deformed_results"]["static"] == []
    assert len(outputs["eos_doc"].volumes) == 11

    fit = _get_fit(outputs["eos_doc"], "birch_murnaghan")
    assert pytest.approx(134.311, abs=1e-2) == fit.B0
    assert pytest.approx(undeformed["atoms"].get_volume(), rel=1e-3) == fit.V0


@pytest.mark.parametrize("pre_relax", [True, False])
def test_eos_flow_run_static(tmp_path, monkeypatch, pre_relax):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu", a=3.5)
    outputs = eos_flow(
        atoms, relax_job, static_job, pre_relax=pre_relax, run_static=True
    )

    # The bulk structure gets a static calculation, after a cell relaxation if
    # pre_relax is True
    undeformed = outputs["undeformed_result"]
    assert undeformed["name"] == "EMT Static"
    assert (
        undeformed["atoms"].get_volume() == pytest.approx(atoms.get_volume())
    ) is not pre_relax

    # With run_static, fixed-cell relaxations then statics on the strained structures
    assert len(outputs["deformed_results"]["relax"]) == 10
    assert len(outputs["deformed_results"]["static"]) == 10
