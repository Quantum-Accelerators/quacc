from __future__ import annotations

from ase import Atoms
from ase.build import bulk, molecule

from quacc import Remove, change_settings
from quacc.recipes.aims import _base, core
from quacc.utils.kpts import kspacing_to_kpts


def _record_calls(result):
    calls = []

    def recorder(*args, **kwargs):
        calls.append((args, kwargs))
        return result

    return recorder, calls


def test_prep_calculator_periodic(tmp_path):
    atoms = bulk("Cu")

    with change_settings(
        {
            "AIMS_BIN": "aims.custom.x",
            "AIMS_PARALLEL_CMD": "mpirun -np 4",
            "AIMS_SPECIES_DEFAULTS": tmp_path,
        }
    ):
        calc = _base.prep_calculator(
            atoms,
            calc_defaults={"kspacing": 0.1, "species_dir": "tight", "xc": "pbe"},
            calc_swaps={"xc": "pbesol"},
        )

    assert calc.parameters["k_grid"] == kspacing_to_kpts(0.1, atoms)
    assert calc.parameters["xc"] == "pbesol"
    assert "kspacing" not in calc.parameters
    assert calc.profile.command == "mpirun -np 4 aims.custom.x"
    assert calc.profile.default_species_directory == str(tmp_path / "tight")


def test_prep_calculator_aperiodic_and_magnetic():
    atoms = molecule("H2")
    atoms.set_initial_magnetic_moments([1.0, 0.0])

    calc = _base.prep_calculator(
        atoms,
        calc_defaults={
            "kspacing": 0.1,
            "k_grid": [2, 2, 2],
            "k_grid_density": 3.0,
            "output_level": "normal",
        },
        calc_swaps={"output_level": Remove},
    )

    assert calc.parameters == {"spin": "collinear"}


def test_prep_calculator_respects_explicit_spin_and_k_grid():
    atoms = bulk("Si")
    atoms.set_initial_magnetic_moments([1.0] * len(atoms))

    calc = _base.prep_calculator(
        atoms, calc_defaults={"kspacing": 0.1, "k_grid": [3, 4, 5], "spin": "none"}
    )

    assert calc.parameters["kspacing"] == 0.1
    assert calc.parameters["k_grid"] == [3, 4, 5]
    assert calc.parameters["spin"] == "none"


def test_run_and_summarize(monkeypatch):
    atoms = molecule("H2")
    final_atoms = atoms.copy()
    calc = object()
    summary = {"energy": -1.0}
    calls = {}

    def prep_calculator(*args, **kwargs):
        calls["prep_calculator"] = (args, kwargs)
        return calc

    class FakeRunner:
        def __init__(self, *args, **kwargs):
            calls["runner_init"] = (args, kwargs)

        def run_calc(self):
            calls["run_calc"] = True
            return final_atoms

    class FakeSummarize:
        def __init__(self, *args, **kwargs):
            calls["summarize_init"] = (args, kwargs)

        def run(self, *args, **kwargs):
            calls["summarize_run"] = (args, kwargs)
            return summary

    monkeypatch.setattr(_base, "prep_calculator", prep_calculator)
    monkeypatch.setattr(_base, "Runner", FakeRunner)
    monkeypatch.setattr(_base, "Summarize", FakeSummarize)

    result = _base.run_and_summarize(
        atoms,
        calc_defaults={"xc": "pbe"},
        calc_swaps={"spin": "none"},
        additional_fields={"name": "test"},
        copy_files="source",
    )

    assert result == summary
    assert calls["prep_calculator"] == (
        (atoms,),
        {"calc_defaults": {"xc": "pbe"}, "calc_swaps": {"spin": "none"}},
    )
    assert calls["runner_init"] == ((atoms, calc), {"copy_files": "source"})
    assert calls["run_calc"] is True
    assert calls["summarize_init"] == (
        (),
        {"move_magmoms": True, "additional_fields": {"name": "test"}},
    )
    assert calls["summarize_run"] == ((final_atoms, atoms), {})


def test_static_job_defaults_and_overrides(monkeypatch):
    atoms = bulk("Cu")
    run_and_summarize, calls = _record_calls({"done": True})
    monkeypatch.setattr(core, "run_and_summarize", run_and_summarize)
    monkeypatch.setattr(core, "check_is_metal", lambda _atoms: True)

    result = core.static_job(
        atoms,
        species_defaults="tight",
        kspacing=0.05,
        spin="none",
        additional_fields={"tag": "custom"},
        copy_files="source",
        xc="pbesol",
    )

    assert result == {"done": True}
    assert len(calls) == 1
    kwargs = calls[0][1]
    assert kwargs["calc_defaults"] == core.BASE_SET_METAL | {
        "species_dir": "tight",
        "kspacing": 0.05,
        "spin": "none",
    }
    assert kwargs["calc_swaps"] == {"xc": "pbesol"}
    assert kwargs["additional_fields"] == {"name": "FHI-aims Static", "tag": "custom"}
    assert kwargs["copy_files"] == "source"


def test_static_job_nonmetal_and_agnostic(monkeypatch):
    atoms = bulk("Si")
    run_and_summarize, calls = _record_calls({})
    monkeypatch.setattr(core, "run_and_summarize", run_and_summarize)
    monkeypatch.setattr(core, "check_is_metal", lambda _atoms: False)

    core.static_job(atoms)
    assert calls[-1][1]["calc_defaults"] == (
        core.BASE_SET_NON_METAL
        | {"species_dir": "intermediate", "kspacing": core.KSPACING_NON_METAL}
    )

    core.static_job(atoms, agnostic_params=True)
    assert calls[-1][1]["calc_defaults"] == (
        core.BASE_SET_AGNOSTIC
        | {"species_dir": "intermediate", "kspacing": core.KSPACING_AGNOSTIC}
    )


def test_relax_job_defaults_and_cell_relaxation(monkeypatch):
    atoms = Atoms("H", cell=[8, 8, 8], pbc=True)
    run_and_summarize, calls = _record_calls({"done": True})
    monkeypatch.setattr(core, "run_and_summarize", run_and_summarize)
    monkeypatch.setattr(core, "check_is_metal", lambda _atoms: False)

    result = core.relax_job(
        atoms,
        relax_cell=True,
        additional_fields={"name": "custom name"},
        force_occupation_type=True,
    )

    assert result == {"done": True}
    assert len(calls) == 1
    kwargs = calls[0][1]
    assert kwargs["calc_defaults"] == core.BASE_SET_NON_METAL | {
        "species_dir": "light",
        "kspacing": core.KSPACING_NON_METAL,
        "relax_geometry": "bfgs 1E-2",
        "relax_unit_cell": "full",
    }
    assert kwargs["calc_swaps"] == {"force_occupation_type": True}
    assert kwargs["additional_fields"] == {"name": "custom name"}


def test_relax_job_agnostic_without_cell_relaxation(monkeypatch):
    atoms = molecule("H2")
    run_and_summarize, calls = _record_calls({})
    monkeypatch.setattr(core, "run_and_summarize", run_and_summarize)

    core.relax_job(atoms, agnostic_params=True, kspacing=0.08, spin="collinear")

    defaults = calls[0][1]["calc_defaults"]
    assert defaults == core.BASE_SET_AGNOSTIC | {
        "species_dir": "light",
        "kspacing": 0.08,
        "spin": "collinear",
        "relax_geometry": "bfgs 1E-2",
    }
    assert "relax_unit_cell" not in defaults
