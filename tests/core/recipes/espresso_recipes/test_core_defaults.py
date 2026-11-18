from __future__ import annotations

from copy import deepcopy

from ase.build import bulk, molecule

from quacc.recipes.espresso import core


def test_jobs_do_not_mutate_base_sets(monkeypatch):
    captured = []

    def capture(atoms, **kwargs):
        captured.append(kwargs["calc_defaults"])
        return {}

    monkeypatch.setattr(core, "run_and_summarize", capture)
    monkeypatch.setattr(core, "run_and_summarize_opt", capture)

    base_set_metal = deepcopy(core.BASE_SET_METAL)
    base_set_non_metal = deepcopy(core.BASE_SET_NON_METAL)

    core.relax_job(bulk("Cu"), relax_cell=True)
    core.ase_relax_job(bulk("Cu"), relax_cell=True)
    core.static_job(bulk("Cu"))
    core.static_job(molecule("H2O"))

    # Each job gets its own defaults, so later (or concurrent) jobs cannot
    # change the settings of an earlier one
    assert captured[0]["input_data"]["control"] == {"calculation": "vc-relax"}
    assert captured[1]["input_data"]["control"]["tstress"] is True
    assert captured[2]["input_data"]["control"] == {"calculation": "scf"}
    assert captured[3]["input_data"]["control"] == {"calculation": "scf"}

    assert base_set_metal == core.BASE_SET_METAL
    assert base_set_non_metal == core.BASE_SET_NON_METAL
