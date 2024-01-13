"""Base jobs for espresso."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from ase.io.espresso import Namelist

from quacc import SETTINGS
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
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
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
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.
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

    bin_path = str(
        Path(SETTINGS.ESPRESSO_BIN_DIR, SETTINGS.ESPRESSO_BINARIES[binary])
    )
    pseudo_path = (
        calc_flags.get("input_data", {})
        .get("control", {})
        .get("pseudo_dir", str(SETTINGS.ESPRESSO_PSEUDO))
    )

    atoms.calc = Espresso(
        input_atoms=atoms,
        preset=preset,
        template=template,
        profile=EspressoProfile(
            bin_path, pseudo_path, argv=SETTINGS.ESPRESSO_PARALLEL_CMD
        ),
        **calc_flags,
    )

    geom_file = template.outputname if template.binary == "pw" else None
    final_atoms = run_calc(atoms, geom_file=geom_file, copy_files=copy_files)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields=additional_fields
    )
