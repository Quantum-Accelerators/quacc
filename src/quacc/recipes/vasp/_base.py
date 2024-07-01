"""Core recipes for VASP."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.vasp import Vasp
from quacc.runners.ase import Runner
from quacc.schemas.vasp import summarize_vasp_opt_run, vasp_summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.runners.ase import OptParams
    from quacc.schemas._aliases.vasp import VaspASEOptSchema, VaspSchema
    from quacc.types import Filenames, SourceDirectory


def run_and_summarize(
    atoms: Atoms,
    preset: str | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> VaspSchema:
    """
    Base job function for VASP recipes.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    calc_defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].
    report_mp_corrections
        Whether to report the Materials Project corrections in the results.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    calc = Vasp(atoms, preset=preset, **calc_flags)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return vasp_summarize_run(
        final_atoms,
        report_mp_corrections=report_mp_corrections,
        additional_fields=additional_fields,
    )


def run_and_summarize_opt(
    atoms: Atoms,
    preset: str | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> VaspASEOptSchema:
    """
    Base job function for VASP recipes with ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    calc_defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].
    opt_defaults
        Default arguments for the ASE optimizer.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.Runner.run_opt][]
    report_mp_corrections
        Whether to report the Materials Project corrections in the results.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    VaspASEOptSchema
        Dictionary of results
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    calc = Vasp(atoms, preset=preset, **calc_flags)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_flags)

    return summarize_vasp_opt_run(
        dyn,
        report_mp_corrections=report_mp_corrections,
        additional_fields=additional_fields,
    )
