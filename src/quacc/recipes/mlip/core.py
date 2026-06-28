"""Core recipes for universal machine-learned interatomic potentials."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.mlip._base import pick_calculator
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.types import OptParams, OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms,
    library: Literal["fairchem", "matcalc", "rootstock"],
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    library
        MLIP library to use:
        - `fairchem` passes `**calc_kwargs` to `FAIRChemCalculator.from_model_checkpoint()`
        - `matcalc` passes `**calc_kwargs` to `matcalc.load_fp()`
        - `rootstock` passes `**calc_kwargs` to `rootstock.RootstockCalculator()`
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the underlying MLIP library.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc = pick_calculator(library, **calc_kwargs)
    final_atoms = Runner(atoms, calc).run_calc()
    return Summarize(
        additional_fields={"name": f"{library} Static"} | (additional_fields or {})
    ).run(final_atoms, atoms)


@job
def relax_job(
    atoms: Atoms,
    library: Literal["fairchem", "matcalc", "rootstock"],
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    library
        MLIP library to use:
        - `fairchem` passes `**calc_kwargs` to `FAIRChemCalculator.from_model_checkpoint()`
        - `matcalc` passes `**calc_kwargs` to `matcalc.load_fp()`
        - `rootstock` passes `**calc_kwargs` to `rootstock.RootstockCalculator()`
    relax_cell
        Whether to relax the cell.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the underlying MLIP library.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    opt_defaults = {"fmax": 0.01}
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    calc = pick_calculator(library, **calc_kwargs)

    dyn = Runner(atoms, calc).run_opt(relax_cell=relax_cell, **opt_flags)

    return Summarize(
        additional_fields={"name": f"{library} Relax"} | (additional_fields or {})
    ).opt(dyn)
