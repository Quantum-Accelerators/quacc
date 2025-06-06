"""Utilities for k-point handling."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import PmgKpts


def convert_pmg_kpts(
    pmg_kpts: PmgKpts, input_atoms: Atoms, force_gamma: bool = False
) -> tuple[list[int], bool]:
    """
    Shortcuts for pymatgen k-point generation schemes.

    Parameters
    ----------
    pmg_kpts
        The pmg_kpts kwargs. Has the following options:

        - {"line_density": float}. This will call
        `pymatgen.symmetry.bandstructure.HighSymmKpath`
            with `path_type="latimer_munro"`. The `line_density` value will be
            set in the `.get_kpoints` attribute.

        - {"kppvol": float}. This will call
        `pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_vol`
            with the given value for `kppvol`.

        - {"kppa": float}. This will call
        `pymatgen.io.vasp.inputs.Kpoints.automatic_density`
            with the given value for `kppa`.

        - {"length_densities": [float, float, float]}. This will call
        `pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_lengths`
            with the given value for `length_densities`.

        If multiple options are specified, the most dense k-point scheme will be
        chosen.
    input_atoms
        The input atoms.
    force_gamma
        Force gamma-centered k-points.

    Returns
    -------
    kpts
        The generated k-points.
    gamma
        Whether the k-points are gamma-centered.
    """
    struct = AseAtomsAdaptor.get_structure(input_atoms)

    if pmg_kpts.get("line_density"):
        kpath = HighSymmKpath(
            struct,
            path_type="latimer_munro",
            has_magmoms=np.any(struct.site_properties.get("magmom", None)),
        )
        kpts, _ = kpath.get_kpoints(
            line_density=pmg_kpts["line_density"], coords_are_cartesian=False
        )
        kpts = np.stack(kpts)
        gamma = False

    else:
        max_pmg_kpts: PmgKpts = None
        for k, v in pmg_kpts.items():
            if k == "kppvol":
                pmg_kpts = Kpoints.automatic_density_by_vol(
                    struct, v, force_gamma=force_gamma
                )
            elif k == "kppa":
                pmg_kpts = Kpoints.automatic_density(struct, v, force_gamma=force_gamma)
            elif k == "length_densities":
                pmg_kpts = Kpoints.automatic_density_by_lengths(
                    struct, v, force_gamma=force_gamma
                )
            else:
                msg = f"Unsupported k-point generation scheme: {pmg_kpts}."
                raise ValueError(msg)

            max_pmg_kpts = (
                pmg_kpts
                if (
                    not max_pmg_kpts
                    or np.prod(pmg_kpts.kpts[0]) >= np.prod(max_pmg_kpts.kpts[0])
                )
                else max_pmg_kpts
            )

        kpts = [int(k) for k in max_pmg_kpts.kpts[0]]
        gamma = max_pmg_kpts.style.name.lower() == "gamma"

    return kpts, gamma


def bandgap_to_kspacing(bandgap: float) -> float:
    """
    Takes the bandgap energy and computes the required KSPACING value.
    Refer to https://drive.google.com/file/d/1fUUx0wrrtMRcSss5yv3NiQuC7J5IiEKL/view

    Parameters
    ----------
    bandgap
        The bandgap of the material in eV.

    Returns
    ----------
    kspacing
        Value of the KSPACING INCAR tag in inverse angstroms.
    """

    deltak_min = 0.2
    deltak_max = 0.45
    a = 0.9
    b = 2.35
    c = 8.0

    delta = a * (bandgap - b)
    return 0.5 * (
        deltak_min
        + deltak_max
        + (deltak_max - deltak_min) * delta / ((1 + delta**c) ** (1 / c))
    )
