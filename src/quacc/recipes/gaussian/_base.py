"""Base jobs for Gaussian."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.gaussian import Gaussian

from quacc import get_settings
from quacc.runners.ase import Runner
from quacc.schemas.cclib import CclibSummarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, SourceDirectory, cclibSchema

_LABEL = "Gaussian"
LOG_FILE = f"{_LABEL}.log"


def run_and_summarize(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> cclibSchema:
    """
    Base job function for carrying out Gaussian recipes.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        Default parameters for the calculator.
    calc_swaps
        Dictionary of custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results, as specified in
        [quacc.schemas.cclib.CclibSummarize.run][]
    """
    settings = get_settings()
    gaussian_cmd = f"{settings.GAUSSIAN_CMD} < {_LABEL}.com > {LOG_FILE}"
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    calc = Gaussian(command=gaussian_cmd, label=_LABEL, **calc_flags)
    atoms = Runner(atoms, calc, copy_files=copy_files).run_calc(geom_file=LOG_FILE)

    return CclibSummarize(LOG_FILE, additional_fields=additional_fields).run(atoms)
