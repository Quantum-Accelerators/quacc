"""Base jobs for espresso."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.atoms import Atoms
from ase.io.espresso import Namelist
from ase.io.espresso_namelist.keys import ALL_KEYS

from quacc.calculators.espresso.espresso import (
    Espresso,
    EspressoProfile,
    EspressoTemplate,
)
from quacc.calculators.espresso.utils import (
    prepare_copy_files,
    remove_conflicting_kpts_kspacing,
)
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from quacc.runners.ase import OptParams
    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


def run_and_summarize(
    atoms: Atoms | None = None,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    parallel_info: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
) -> RunSchema:
    """
    Base function to carry out espresso recipes.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Name of the preset to use
    template
        EspressoTemplate to use
    profile
        EspressoProfile to use
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.espresso.Espresso][] calculator.
    parallel_info
        Dictionary of parallelization information.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms = prepare_atoms(
        atoms=atoms,
        preset=preset,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        calc_swaps=calc_swaps,
        parallel_info=parallel_info,
    )

    updated_copy_files = prepare_copy(
        copy_files=copy_files,
        calc_params=atoms.calc.user_calc_params,
        binary=atoms.calc.template.binary,
    )

    geom_file = template.outputname if template.binary == "pw" else None

    final_atoms = run_calc(atoms, geom_file=geom_file, copy_files=updated_copy_files)

    return summarize_run(
        final_atoms, atoms, move_magmoms=True, additional_fields=additional_fields
    )


def run_and_summarize_opt(
    atoms: Atoms | None = None,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    parallel_info: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
) -> RunSchema:
    """
    Base function to carry out espresso recipes with ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Name of the preset to use
    template
        EspressoTemplate to use
    profile
        EspressoProfile to use
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.espresso.Espresso][] calculator.
    opt_defaults
        The default optimization parameters.
    opt_params
        Dictionary of parameters to pass to the optimizer. pass "optimizer"
        to change the optimizer being used. "fmax" and "max_steps" are commonly
        used keywords. See the ASE documentation for more information.
    parallel_info
        Dictionary of parallelization information.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms = prepare_atoms(
        atoms=atoms,
        preset=preset,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        calc_swaps=calc_swaps,
        parallel_info=parallel_info,
    )

    updated_copy_files = prepare_copy(
        copy_files=copy_files,
        calc_params=atoms.calc.user_calc_params,
        binary=atoms.calc.template.binary,
    )

    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    dyn = run_opt(atoms, copy_files=updated_copy_files, **opt_flags)

    return summarize_opt_run(
        dyn, move_magmoms=True, additional_fields=additional_fields
    )


def prepare_atoms(
    atoms: Atoms | None = None,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    parallel_info: dict[str, Any] | None = None,
) -> Atoms:
    """
    Commonly used preparation function to merge parameters
    and attach an Espresso calculator accordingly.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Name of the preset to use
    template
        EspressoTemplate to use
    profile
        EspressoProfile to use
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.espresso.Espresso][] calculator.
    parallel_info
        Dictionary of parallelization information.

    Returns
    -------
    Atoms
        Atoms object with attached Espresso calculator.
    """
    atoms = Atoms() if atoms is None else atoms
    calc_defaults = calc_defaults or {}
    calc_swaps = calc_swaps or {}

    calc_defaults["input_data"] = Namelist(calc_defaults.get("input_data"))
    calc_swaps["input_data"] = Namelist(calc_swaps.get("input_data"))

    binary = template.binary if template else "pw"

    if binary in ALL_KEYS:
        calc_defaults["input_data"].to_nested(binary=binary, **calc_defaults)
        calc_swaps["input_data"].to_nested(binary=binary, **calc_swaps)

    calc_defaults = remove_conflicting_kpts_kspacing(calc_defaults, calc_swaps)

    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    atoms.calc = Espresso(
        input_atoms=atoms,
        preset=preset,
        parallel_info=parallel_info,
        template=template,
        profile=profile,
        **calc_flags,
    )

    return atoms


def prepare_copy(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    calc_params: dict[str, Any] | None = None,
    binary: str = "pw",
) -> dict[SourceDirectory, Filenames] | None:
    """
    Function that will prepare the files to copy.

    Parameters
    ----------
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    calc_params
        The calculator parameters.
    binary
        The binary to use.

    Returns
    -------
    dict
        Dictionary of files to copy.
    """
    if isinstance(copy_files, (str, Path)):
        copy_files = [copy_files]

    if isinstance(copy_files, list):
        exact_files_to_copy = prepare_copy_files(calc_params, binary=binary)
        return {source: exact_files_to_copy for source in copy_files}

    return copy_files
