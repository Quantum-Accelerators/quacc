"""Base jobs for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms
from ase.io.espresso import Namelist

from quacc.calculators.espresso.espresso import (
    Espresso,
    EspressoProfile,
    EspressoTemplate,
)
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from quacc.schemas._aliases.ase import RunSchema


def base_fn(
    atoms: Atoms = None,
    preset: str | None = None,
    template: EspressoTemplate | None = None,
    profile: EspressoProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    parallel_info: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
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
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.
    parallel_info
        Dictionary of parallelization information.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    atoms = Atoms() if atoms is None else atoms

    calc_defaults["input_data"] = Namelist(calc_defaults.get("input_data"))
    calc_swaps["input_data"] = Namelist(calc_swaps.get("input_data"))

    binary = template.binary if template else "pw"

    calc_defaults["input_data"].to_nested(binary=binary, **calc_defaults)
    calc_swaps["input_data"].to_nested(binary=binary, **calc_swaps)

    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    atoms.calc = Espresso(
        input_atoms=atoms,
        preset=preset,
        parallel_info=parallel_info,
        template=template,
        profile=profile,
        **calc_flags,
    )

    geom_file = template.outputname if template.binary == "pw" else None

    final_atoms = run_calc(atoms, geom_file=geom_file, copy_files=copy_files)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields=additional_fields
    )
