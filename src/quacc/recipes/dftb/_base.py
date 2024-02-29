"""Base jobs for DFTB+"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.dftb import Dftb

from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory

LOG_FILE = "dftb.out"
GEOM_FILE = "geo_end.gen"


def base_fn(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> RunSchema:
    """
    Base job function for DFTB+ recipes.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator that would override the
        calculator defaults. Set a value to `quacc.Remove` to remove a pre-existing key
        entirely. For a list of available keys, refer to the
        [ase.calculators.dftb.Dftb][] calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """

    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    atoms.calc = Dftb(**calc_flags)
    final_atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(final_atoms, atoms, additional_fields=additional_fields)
