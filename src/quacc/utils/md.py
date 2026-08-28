"""Utility functions for molecular dynamics."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.md.andersen import Andersen
from ase.md.bussi import Bussi
from ase.md.langevin import Langevin
from ase.md.nose_hoover_chain import MTKNPT, IsotropicMTKNPT, NoseHooverChainNVT
from ase.md.npt import NPT
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen, NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.units import bar, fs

from quacc.utils.dicts import Remove

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.md.md import MolecularDynamics

    from quacc.types import MDEnsemble


def resolve_md_ensemble(
    ensemble: MDEnsemble,
    timestep: float,
    temperature_K: float | None,
    pressure_bar: float | None,
) -> tuple[type[MolecularDynamics], dict[str, Any]]:
    """
    Resolve an ensemble name to an ASE dynamics class and its default kwargs.

    Parameters
    ----------
    ensemble
        Name of the ensemble (case-insensitive). See [quacc.types.MDEnsemble][]
        for the supported names.
    timestep
        Time step in ASE units.
    temperature_K
        Temperature in K, if applicable for the given ensemble.
    pressure_bar
        Pressure in bar, if applicable for the given ensemble.

    Returns
    -------
    tuple[type[MolecularDynamics], dict]
        The ASE dynamics class and its default `dynamics_kwargs`.
    """
    taut = 100 * timestep
    taup = 1000 * timestep
    temperature = {
        "temperature_K": temperature_K if temperature_K is not None else Remove
    }
    pressure = (pressure_bar if pressure_bar is not None else 0.0) * bar

    presets: dict[str, tuple[type[MolecularDynamics], dict[str, Any]]] = {
        "nve": (VelocityVerlet, {}),
        "nvt": (NoseHooverChainNVT, temperature | {"tdamp": taut}),
        "nvt_berendsen": (NVTBerendsen, temperature | {"taut": taut}),
        "nvt_langevin": (Langevin, temperature | {"friction": 1.0e-3}),
        "nvt_andersen": (Andersen, temperature | {"andersen_prob": 1.0e-2}),
        "nvt_bussi": (Bussi, temperature | {"taut": taut}),
        "npt": (
            NPT,
            temperature
            | {"externalstress": pressure, "ttime": 25 * fs, "pfactor": 75.0**2 * fs},
        ),
        "npt_berendsen": (
            NPTBerendsen,
            temperature | {"pressure_au": pressure, "taut": taut, "taup": taup},
        ),
        "npt_inhomogeneous": (
            Inhomogeneous_NPTBerendsen,
            temperature | {"pressure_au": pressure, "taut": taut, "taup": taup},
        ),
        "npt_mtk": (
            MTKNPT,
            temperature | {"pressure_au": pressure, "tdamp": taut, "pdamp": taup},
        ),
        "npt_isotropic_mtk": (
            IsotropicMTKNPT,
            temperature | {"pressure_au": pressure, "tdamp": taut, "pdamp": taup},
        ),
    }

    if ensemble.lower() not in presets:
        msg = (
            f"Unsupported ensemble: {ensemble}. "
            f"Valid options are: {', '.join(presets)}."
        )
        raise ValueError(msg)

    return presets[ensemble.lower()]


def upper_triangular_cell(atoms: Atoms) -> Atoms:
    """
    Rotate the cell to the upper-triangular form required by the ASE `NPT`
    class, using `ase.cell.Cell.standard_form`. The positions and any momenta are rotated
    along with the cell, so the structure is unchanged.

    Parameters
    ----------
    atoms
        Atoms object.

    Returns
    -------
    Atoms
        A copy of the Atoms object with an upper-triangular cell, or the
        original Atoms object if the cell is already upper-triangular.
    """
    if atoms.cell[1, 0] == atoms.cell[2, 0] == atoms.cell[2, 1] == 0.0:
        return atoms

    atoms = atoms.copy()
    cell, rotation = atoms.cell.standard_form(form="upper")
    atoms.set_cell(cell)
    atoms.set_positions(atoms.get_positions() @ rotation.T)
    if atoms.has("momenta"):
        atoms.set_momenta(atoms.get_momenta() @ rotation.T)
    return atoms
