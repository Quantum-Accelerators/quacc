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
    atoms: Atoms, copy_files: list[str] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a basic SCF calculation with ONETEP.

    Parameters
    ----------
    atoms
        The Atoms object.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.onetep.Onetep` calculator.

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
    copy_files: list[str] | None = None,
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
        List of files to copy to the calculation directory. Useful for copying
        files from a previous calculation. This parameter can either be a string
        or a list of strings.

        If a string is provided, it is assumed to be a path to a directory,
        all of the child tree structure of that directory is going to be copied to the
        scratch of this calculation.

        If a list of strings is provided, each string point to a specific file. In this case
        it is important to note that no directory structure is going to be copied, everything
        is copied at the root of the temporary directory.
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

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": LBFGS}

    return base_opt_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "ONETEP ASE Relax"},
        copy_files=copy_files,
    )
