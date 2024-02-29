"""Core recipes for Onetep."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import LBFGS

from quacc import job
from quacc.recipes.onetep._base import base_fn, base_opt_fn
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory

BASE_SET = {
    "keywords": {
        "output_detail": "verbose",
        "do_properties": True,
        "cutoff_energy": "600 eV",
        "task": "SinglePoint",
    }
}


@job
def static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic SCF calculation with ONETEP.

    Parameters
    ----------
    atoms
        The Atoms object.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.onetep.Onetep][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = BASE_SET

    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "ONETEP Static"},
        copy_files=copy_files,
    )


@job
def ase_relax_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    opt_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a structure relaxation with ONETEP using ASE
    external optimizers.

    Parameters
    ----------
    atoms
        The Atoms object.
    opt_params
        Dictionary of parameters to pass to the optimizer. pass "optimizer"
        to change the optimizer being used. "fmax" and "max_steps" are commonly
        used keywords. See the ASE documentation for more information.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Additional keyword arguments to pass to the ONETEP calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = recursive_dict_merge(
        BASE_SET,
        {"keywords": {"write_forces": True, "forces_output_detail": "verbose"}},
    )

    opt_defaults = {"optimizer": LBFGS}

    return base_opt_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "ONETEP ASE Relax"},
        copy_files=copy_files,
    )
