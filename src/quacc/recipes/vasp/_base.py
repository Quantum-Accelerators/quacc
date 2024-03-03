"""Core recipes for VASP."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.vasp import Vasp
from quacc.runners.ase import run_calc
from quacc.schemas.vasp import vasp_summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import VaspSchema
    from quacc.utils.files import Filenames, SourceDirectory


def base_fn(
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
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    atoms.calc = Vasp(atoms, preset=preset, **calc_flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return vasp_summarize_run(
        final_atoms,
        report_mp_corrections=report_mp_corrections,
        additional_fields=additional_fields,
    )
