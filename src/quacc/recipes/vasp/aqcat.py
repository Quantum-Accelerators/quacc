"""AQCat25-compatible slab recipes"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from quacc import job
from quacc.recipes.vasp._base import run_and_summarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import CopyFiles, VaspSchema


@job
def aqcat25_static_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Carry out a static calculation with AQCat25 slab settings.

    By default, use the slab k-point sampling from the AQCat25 paper's
    fairchem implementation: max(1, round(40 / L)) in each in-plane direction,
    where L is the largest absolute Cartesian component of that cell vector
    in angstroms, and one k-point along the third cell vector (vacuum).
    Override this sampling with `kpts` or `kspacing` in `calc_kwargs`.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    # Initial magnetic moments of the spin-polarized elements (AQCat25 paper, Table 6)
    magmoms = {
        "V": 5.0,
        "Cr": 5.0,
        "Mn": 5.0,
        "Fe": 5.0,
        "Co": 5.0,
        "Ni": 5.0,
        "Cu": 1.73,
        "Mo": 5.0,
        "Ru": 2.2,
        "W": 5.0,
        "Os": 2.2,
        "Ce": 5.0,
    }
    spin_polarized = not set(atoms.get_chemical_symbols()).isdisjoint(magmoms)

    calc_defaults = {
        "ibrion": 2,
        "nsw": 0,
        "isif": 0,
        "ispin": 2 if spin_polarized else 1,
        "isym": 0,
        "algo": "normal",
        "ismear": 0,
        "sigma": 0.1,
        "ediffg": -0.03,
        "encut": 500.0,
        "prec": "accurate",
        "potim": 0.5,
        "nelm": 250,
        "lwave": False,
        "lvhar": False,
        "lcharg": False,
        "laechg": False,
        "xc": "rpbe",
        "lasph": False,
        "ediff": 1e-4,
        "symprec": 1e-10,
        "lreal": "auto",
        "pp_version": "54",
        "incar_copilot_mode": "critical",
        "use_custodian": False,
    }
    if spin_polarized:
        calc_defaults |= {"elemental_magmoms": magmoms, "preset_mag_default": 0.0}
    if "kpts" not in calc_kwargs and "kspacing" not in calc_kwargs:
        # Match fairchem's calculate_surface_k_points used by AQCat25.
        cell = atoms.get_cell()
        calc_defaults["kpts"] = (
            max(1, round(40 / np.linalg.norm(cell[0], ord=np.inf))),
            max(1, round(40 / np.linalg.norm(cell[1], ord=np.inf))),
            1,
        )

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "AQCat25 Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )
