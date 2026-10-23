from __future__ import annotations

from pathlib import Path

import pytest
from ase import Atoms
from ase.build import bulk, molecule
from ase.data import atomic_numbers

from quacc import Remove, change_settings
from quacc.recipes.aims import _base, mad


@pytest.fixture
def write_inputs(monkeypatch, tmp_path):
    """Exercise ASE's actual input writer without running a licensed executable."""
    species_dir = tmp_path / "tight"
    species_dir.mkdir()
    captured = {}

    def run(atoms, calc_defaults, calc_swaps, additional_fields, copy_files):
        calc = _base.prep_calculator(atoms, calc_defaults, calc_swaps)
        calc.template.write_input(
            calc.profile, tmp_path, atoms, calc.parameters, ["energy", "forces"]
        )
        captured.update(
            parameters=calc.parameters,
            control=(tmp_path / "control.in").read_text(),
            fields=additional_fields,
            copy_files=copy_files,
            species_dir=Path(calc.profile.default_species_directory),
        )
        return captured

    monkeypatch.setattr(mad, "run_and_summarize", run)
    with change_settings({"AIMS_SPECIES_DEFAULTS": tmp_path}):
        yield species_dir, captured


def _write_species(directory, symbols):
    # Synthetic species content, not redistributed FHI-aims basis sets.
    for symbol in symbols:
        (directory / f"{atomic_numbers[symbol]:02}_{symbol}_default").write_text(
            f"species {symbol}\n"
            "  confined 5 d auto\n"
            "  confined 6 d auto\n"
            "  confined 5 f auto\n"
            "  hydro 4 f 7.4\n"
            "# confined 5 d 1.0\n"
        )


def _control_flags(text):
    flags = []
    for line in text.splitlines():
        tokens = line.split("#", 1)[0].split()
        if tokens and tokens[0] == "species":
            break
        if tokens:
            flags.append(tokens)
    return flags


def test_published_control(write_inputs):
    species_dir, captured = write_inputs
    _write_species(species_dir, ["Cu"])
    mad.static_job(bulk("Cu"))

    expected = sorted(
        _control_flags((Path(__file__).parent / "mad_control.in").read_text()),
        key=lambda row: row[0],
    )
    actual = sorted(_control_flags(captured["control"]), key=lambda row: row[0])
    assert [row[0] for row in actual] == [row[0] for row in expected]
    for actual_row, expected_row in zip(actual, expected, strict=True):
        assert len(actual_row) == len(expected_row)
        for actual_value, expected_value in zip(actual_row, expected_row, strict=True):
            try:
                assert float(actual_value.replace("d", "e")) == float(
                    expected_value.replace("d", "e")
                )
            except ValueError:
                assert actual_value == expected_value
    assert captured["fields"] == {"name": "MAD-1.5/1.6 Static"}
    assert captured["species_dir"] == species_dir


@pytest.mark.parametrize("cell", [None, [10, 10, 10]])
@pytest.mark.parametrize("use_dipole_correction", [True, False])
def test_molecule(write_inputs, cell, use_dipole_correction):
    species_dir, captured = write_inputs
    _write_species(species_dir, ["H"])
    atoms = molecule("H2")
    if cell:
        atoms.set_cell(cell)
    atoms.set_initial_magnetic_moments([1, 1])
    mad.static_job(atoms, use_dipole_correction=use_dipole_correction)
    assert captured["parameters"]["spin"] == "none"
    for key in ("k_grid_density", "compute_analytical_stress", "use_dipole_correction"):
        assert key not in captured["parameters"]


@pytest.mark.parametrize("use_dipole_correction", [True, False])
def test_dipole_correction(write_inputs, use_dipole_correction):
    species_dir, captured = write_inputs
    _write_species(species_dir, ["Si"])
    atoms = Atoms("Si", cell=[5, 5, 20], pbc=[True, True, False])
    mad.static_job(atoms, use_dipole_correction=use_dipole_correction)
    assert (
        captured["parameters"].get("use_dipole_correction", False)
        == use_dipole_correction
    )
    assert captured["parameters"]["k_grid_density"] == 8
    assert captured["parameters"]["compute_analytical_stress"] is True


@pytest.mark.parametrize(
    ("symbol", "shell"),
    [
        ("La", None),
        ("Ce", None),
        ("Pr", "5"),
        ("Nd", "5"),
        ("Yb", "5"),
        ("Lu", None),
        ("Np", None),
        ("Pu", "6"),
        ("Am", "6"),
        ("No", "6"),
    ],
)
def test_confined_basis(write_inputs, symbol, shell):
    species_dir, captured = write_inputs
    _write_species(species_dir, [symbol, "H"])
    originals = {file.name: file.read_text() for file in species_dir.iterdir()}
    mad.static_job(Atoms([symbol, "H"], positions=[[0, 0, 0], [0, 0, 2]]))
    control = captured["control"]
    species_text = control.split(f"species {symbol}\n", 1)[1].split("species ", 1)[0]
    for orbital in ("5", "6"):
        assert (f"  confined {orbital} d auto\n" in species_text) == (orbital != shell)
    assert "  confined 5 f auto\n" in species_text
    assert "  hydro 4 f 7.4\n" in species_text
    assert "# confined 5 d 1.0\n" in species_text
    assert "species H\n" in control
    assert originals == {file.name: file.read_text() for file in species_dir.iterdir()}
    if shell:
        assert not captured["species_dir"].exists()


@pytest.mark.parametrize(
    "kpoints", [{"k_grid": [2, 3, 4]}, {"kpts": [2, 3, 4]}, {"kspacing": 0.2}]
)
def test_overrides(write_inputs, kpoints):
    species_dir, captured = write_inputs
    _write_species(species_dir, ["Cu"])
    atoms = bulk("Cu")
    mad.static_job(
        atoms,
        xc="pbe",
        use_dipole_correction=False,
        override_warning_libxc=Remove,
        additional_fields={"name": "custom", "tag": 1},
        copy_files=[{"source": "source", "filenames": "*"}],
        **kpoints,
    )
    assert "k_grid_density" not in captured["parameters"]
    assert "override_warning_libxc" not in captured["parameters"]
    assert captured["parameters"]["xc"] == "pbe"
    assert "use_dipole_correction" not in captured["parameters"]
    assert captured["fields"] == {"name": "custom", "tag": 1}
    assert captured["copy_files"] == [{"source": "source", "filenames": "*"}]


def test_custom_species_path(write_inputs, tmp_path):
    _, captured = write_inputs
    custom_dir = tmp_path / "custom"
    custom_dir.mkdir()
    _write_species(custom_dir, ["Pr"])
    mad.static_job(Atoms("Pr"), species_dir=custom_dir)
    assert "  confined 5 d auto\n" not in captured["control"]
    assert "  confined 6 d auto\n" in captured["control"]


def test_removing_absent_grid_preserves_density(write_inputs):
    species_dir, captured = write_inputs
    _write_species(species_dir, ["Cu"])
    mad.static_job(bulk("Cu"), k_grid=Remove)
    assert captured["parameters"]["k_grid_density"] == 8


def test_species_cleanup_on_failure(tmp_path):
    _write_species(tmp_path, ["Pr"])
    directories = []

    def fail():
        with mad._species_directory(Atoms("Pr"), tmp_path) as directory:
            directories.append(directory)
            assert directory.exists()
            raise RuntimeError("failed calculation")

    with pytest.raises(RuntimeError, match="failed calculation"):
        fail()
    directory = directories[0]
    assert not directory.exists()
