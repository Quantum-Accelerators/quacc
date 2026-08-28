from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk
from ase.md.andersen import Andersen
from ase.md.bussi import Bussi
from ase.md.langevin import Langevin
from ase.md.nose_hoover_chain import MTKNPT, IsotropicMTKNPT, NoseHooverChainNVT
from ase.md.npt import NPT
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen, NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.units import bar, fs

from quacc import Remove
from quacc.utils.md import resolve_md_ensemble, upper_triangular_cell


@pytest.mark.parametrize(
    ("ensemble", "expected_class"),
    [
        ("nve", VelocityVerlet),
        ("nvt", NoseHooverChainNVT),
        ("nvt_berendsen", NVTBerendsen),
        ("nvt_langevin", Langevin),
        ("nvt_andersen", Andersen),
        ("nvt_bussi", Bussi),
        ("npt", NPT),
        ("npt_berendsen", NPTBerendsen),
        ("npt_inhomogeneous", Inhomogeneous_NPTBerendsen),
        ("npt_mtk", MTKNPT),
        ("npt_isotropic_mtk", IsotropicMTKNPT),
    ],
)
def test_resolve_md_ensemble_classes(ensemble, expected_class):
    dynamics, _ = resolve_md_ensemble(ensemble, 1.0, 300, None)
    assert dynamics is expected_class


def test_resolve_md_ensemble_is_case_insensitive():
    dynamics, _ = resolve_md_ensemble("NVT_Langevin", 1.0, 300, None)
    assert dynamics is Langevin


def test_resolve_md_ensemble_coupling_times_scale_with_timestep():
    _, kwargs = resolve_md_ensemble("npt_berendsen", 2.0, 300, 1.0)
    assert kwargs["taut"] == pytest.approx(200 * fs)
    assert kwargs["taup"] == pytest.approx(2000 * fs)

    _, kwargs = resolve_md_ensemble("nvt", 2.0, 300, None)
    assert kwargs["tdamp"] == pytest.approx(200 * fs)


def test_resolve_md_ensemble_pressure_routing():
    # ase.md.npt.NPT takes `externalstress`; the Berendsen and MTK families
    # take `pressure_au`.
    _, kwargs = resolve_md_ensemble("npt", 1.0, 300, 2.0)
    assert kwargs["externalstress"] == pytest.approx(2.0 * bar)
    assert "pressure_au" not in kwargs

    _, kwargs = resolve_md_ensemble("npt_berendsen", 1.0, 300, 2.0)
    assert kwargs["pressure_au"] == pytest.approx(2.0 * bar)
    assert "externalstress" not in kwargs

    # The pressure defaults to 0 bar for the NPT-family ensembles.
    _, kwargs = resolve_md_ensemble("npt_mtk", 1.0, 300, None)
    assert kwargs["pressure_au"] == 0.0


def test_resolve_md_ensemble_temperature_handling():
    # An explicit 0 K is honored; only None means unset.
    _, kwargs = resolve_md_ensemble("nvt_langevin", 1.0, 0, None)
    assert kwargs["temperature_K"] == 0

    _, kwargs = resolve_md_ensemble("nvt_langevin", 1.0, None, None)
    assert kwargs["temperature_K"] is Remove

    # NVE ignores the temperature and pressure entirely.
    _, kwargs = resolve_md_ensemble("nve", 1.0, 300, 1.0)
    assert kwargs == {}


def test_resolve_md_ensemble_unknown():
    with pytest.raises(ValueError, match="Unsupported ensemble"):
        resolve_md_ensemble("nvt_gibberish", 1.0, 300, None)


def test_upper_triangular_cell():
    # The primitive fcc cell is not upper-triangular.
    atoms = bulk("Cu")
    old_cell = atoms.cell.copy()

    transformed = upper_triangular_cell(atoms)
    new_cell = transformed.cell

    assert transformed is not atoms
    assert new_cell[1, 0] == new_cell[2, 0] == new_cell[2, 1] == 0.0
    assert new_cell.cellpar() == pytest.approx(old_cell.cellpar())
    # The input Atoms object should not be mutated.
    assert atoms.cell == pytest.approx(old_cell[:])

    # An already upper-triangular cell is returned as-is.
    cubic = bulk("Cu", cubic=True)
    assert upper_triangular_cell(cubic) is cubic


def test_upper_triangular_cell_rotates_momenta():
    atoms = bulk("Cu") * (2, 2, 2)
    rng = np.random.default_rng(seed=42)
    atoms.set_momenta(rng.normal(size=(len(atoms), 3)))
    speeds = np.linalg.norm(atoms.get_momenta(), axis=1)

    transformed = upper_triangular_cell(atoms)

    # The momenta are rotated along with the cell: the magnitudes are
    # unchanged, and the scaled (fractional) momenta are preserved.
    assert np.linalg.norm(transformed.get_momenta(), axis=1) == pytest.approx(speeds)
    assert transformed.get_momenta() @ np.linalg.inv(transformed.cell) == pytest.approx(
        atoms.get_momenta() @ np.linalg.inv(atoms.cell)
    )
