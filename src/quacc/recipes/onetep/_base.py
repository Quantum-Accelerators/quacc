"""Base jobs for Onetep."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.onetep import Onetep, OnetepProfile

from quacc import SETTINGS
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema


def base_fn(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
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
        keys, refer to the `ase.calculators.onetep.Onetep` calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    atoms.calc = Onetep(
        calc_defaults=calc_defaults,
        pseudo_path=str(SETTINGS.ONETEP_PP_PATH) if SETTINGS.ONETEP_PP_PATH else ".",
        parallel_info=SETTINGS.ONETEP_PARALLEL_CMD,
        profile=OnetepProfile(
            SETTINGS.ONETEP_CMD
        ),  # TODO: If the ASE merge is successful, we need to change ONETEP_PARALLEL_CMD to a list[str] and remove parallel info.
        # If we also have access to post_args we can point not to the binary but to the launcher which takes -t nthreads as a post_args
        **calc_flags,
    )

    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(final_atoms, atoms, additional_fields=additional_fields)
