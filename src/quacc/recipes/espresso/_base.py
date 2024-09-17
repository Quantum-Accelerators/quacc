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
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from quacc.types import Filenames, OptParams, RunSchema, SourceDirectory


def run_and_summarize(
    atoms: Atoms | None = None,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
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
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    atoms = Atoms() if atoms is None else atoms
    calc = prepare_calc(
        atoms,
        preset=preset,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        calc_swaps=calc_swaps,
    )

    updated_copy_files = prepare_copy(
        copy_files=copy_files,
        calc_params=calc.user_calc_params,
        binary=calc.template.binary,
    )

    geom_file = template.outputname if template and template.binary == "pw" else None

    final_atoms = Runner(atoms, calc, copy_files=updated_copy_files).run_calc(
        geom_file=geom_file
    )

    return Summarize(move_magmoms=True, additional_fields=additional_fields).run(
        final_atoms, atoms
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
    atoms = Atoms() if atoms is None else atoms
    calc = prepare_calc(
        atoms,
        preset=preset,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        calc_swaps=calc_swaps,
    )

    updated_copy_files = prepare_copy(
        copy_files=copy_files,
        calc_params=calc.user_calc_params,
        binary=calc.template.binary,
    )

    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    dyn = Runner(atoms, calc, copy_files=updated_copy_files).run_opt(**opt_flags)

    return Summarize(move_magmoms=True, additional_fields=additional_fields).opt(dyn)


def prepare_calc(
    atoms: Atoms,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
) -> Espresso:
    """
    Commonly used preparation function to merge parameters
    and create an Espresso calculator accordingly.

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

    Returns
    -------
    Espresso
        The Espresso calculator.
    """
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

    return Espresso(
        input_atoms=atoms,
        preset=preset,
        template=template,
        profile=profile,
        **calc_flags,
    )


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
    if isinstance(copy_files, str | Path):
        copy_files = [copy_files]

    if isinstance(copy_files, list):
        exact_files_to_copy = prepare_copy_files(calc_params, binary=binary)
        return {source: exact_files_to_copy for source in copy_files}

    return copy_files
