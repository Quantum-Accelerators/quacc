"""Core recipes for the tblite code"""
from __future__ import annotations

from typing import TYPE_CHECKING

import covalent as ct
from ase.optimize import FIRE
from monty.dev import requires

from quacc.schemas.ase import summarize_opt_run
from quacc.schemas.atoms import fetch_atoms
from quacc.util.calc import run_ase_opt

try:
    import matgl
    from matgl.ext.ase import M3GNetCalculator

except ImportError:
    matgl = None

if TYPE_CHECKING:
    from ase import Atoms
    from torch import Tensor

    from quacc.schemas.ase import OptSchema

@ct.electron
@requires(matgl, "matgl must be installed. Try pip install matgl")
def relax_job(
    atoms: Atoms | dict,
    model: str = "M3GNet-MP-2021.2.8-PES",
    state_attr: Tensor | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    model
        Pre-trained M3GNet model to use for the calculation
    state_attr
        State attribute to use for the calculation
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """
    atoms = fetch_atoms(atom)
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = opt_defaults | opt_swaps

    potential = matgl.load_model(model)
    atoms.calc = M3GNetCalculator(potential=potential, state_attr=state_attr)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "M3GNet Relax"})
