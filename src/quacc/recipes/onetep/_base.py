"""Base jobs for Onetep."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.onetep import Onetep, OnetepProfile

from quacc import get_settings
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, RunSchema, SourceDirectory


def run_and_summarize(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> RunSchema:
    """
    Base function to carry out Onetep recipes.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.onetep.Onetep][] calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    calc = prep_calculator(calc_defaults=calc_defaults, calc_swaps=calc_swaps)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return Summarize(additional_fields=additional_fields).run(final_atoms, atoms)


def run_and_summarize_opt(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> RunSchema:
    """
    Base function to carry out Onetep recipes with ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.onetep.Onetep][] calculator.
    opt_defaults
        The default optimization parameters.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)
    calc = prep_calculator(calc_defaults=calc_defaults, calc_swaps=calc_swaps)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_flags)

    return Summarize(additional_fields=additional_fields).opt(dyn)


def prep_calculator(
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
) -> Onetep:
    """
    Prepare the Onetep calculator.

    Parameters
    ----------
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.onetep.Onetep][] calculator.

    Returns
    -------
    Onetep
        The Onetep calculator.
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
    settings = get_settings()

    return Onetep(
        profile=OnetepProfile(f"{settings.ONETEP_CMD}", str(settings.ONETEP_PP_PATH)),
        **calc_flags,
    )
