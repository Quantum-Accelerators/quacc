"""
Core recipes for the CHGNet code
"""
from __future__ import annotations

from copy import deepcopy

import covalent as ct
from ase.atoms import Atoms
from ase.constraints import ExpCellFilter
from monty.dev import requires

from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc

try:
    from chgnet.model.dynamics import CHGNetCalculator as CHGNet
    from chgnet.model.model import CHGNet as model
except ImportError:
    CHGNet = None


@ct.electron
@requires(CHGNet, "CHGNet must be installed. Try pip install chgnet")
def static_job(
    atoms: Atoms,
    model_path: str | None = None,
    chgnet_kwargs: dict | None = None,
    model_kwargs: dict | None = None,
) -> dict:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    model_path
        Path to the model file
    chgnet_kwargs
        Dictionary of custom kwargs for the CHGNet calculator
    model_kwargs
        Dictionary of custom kwargs for the CHGNet model

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """
    chgnet_kwargs = chgnet_kwargs or {}
    model_kwargs = model_kwargs or {}
    input_atoms = deepcopy(atoms)

    if model_path:
        model.from_file(model_path, **model_kwargs)
    else:
        model.load()
    atoms.calc = CHGNet(model=model, **chgnet_kwargs)
    atoms = run_calc(atoms)
    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "CHGNet Static"},
    )


@ct.electron
@requires(CHGNet, "CHGNet must be installed. Try pip install chgnet")
def relax_job(
    atoms: Atoms,
    relax_cell: bool = True,
    model_path: str | None = None,
    chgnet_kwargs: dict | None = None,
    model_kwargs: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    model_path
        Path to the model file
    chgnet_kwargs
        Dictionary of custom kwargs for the CHGNet calculator
    model_kwargs
        Dictionary of custom kwargs for the CHGNet model
    opt_swaps
        Dictionary of custom kwargs for the optimizer

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """
    chgnet_kwargs = chgnet_kwargs or {}
    model_kwargs = model_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "FIRE"}
    opt_flags = opt_defaults | opt_swaps

    if relax_cell:
        atoms = ExpCellFilter(atoms)

    if model_path:
        model.from_file(model_path, **model_kwargs)
    else:
        model.load()
    atoms.calc = CHGNet(model=model, **chgnet_kwargs)
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "CHGNet Relax"})
