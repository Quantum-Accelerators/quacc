"""Core recipes for Onetep."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import LBFGS

from quacc import job
from quacc.recipes.onetep._base import run_and_summarize, run_and_summarize_opt
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.runners.ase import OptParams
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

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "ONETEP Static"},
        copy_files=copy_files,
    )


@job
def ase_relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a structure relaxation with ONETEP using ASE
    external optimizers.

    Parameters
    ----------
    atoms
        The Atoms object.
    relax_cell
        True if a volume relaxation should be performed. False if only the positions
        should be updated.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.run_opt][].
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

    opt_defaults = {"optimizer": LBFGS, "relax_cell": relax_cell}

    return run_and_summarize_opt(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "ONETEP ASE Relax"},
        copy_files=copy_files,
    )
