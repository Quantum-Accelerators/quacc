"""Core recipes for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.espresso import (Espresso,
                                      EspressoProfile)

from quacc import SETTINGS, job
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import merge_dicts

ESPRESSO_CMD = f"{SETTINGS.ESPRESSO_CMD}"

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    copy_files: list[str] | None = None,
    **kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **kwargs
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "ecutwfc": 40,
                "ecutrho": 160,
                "nspin": 1,
            }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        "ecutwfc": 40,
        "ecutrho": 160,
        "nspin": 1,
    }
    return _base_job(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=kwargs,
        additional_fields={"name": "espresso Static"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms,
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
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
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

    # No reason to type pseudopotentials=pseudopotentials everywhere
    # if something like this is possible?
    pseudopotentials = atoms.info.get("pseudopotentials", None)
    # No reason to add a env variable, pseudo path is already
    # managed by $ESPRESSO_PSEUDO.

    flags = merge_dicts(calc_defaults, calc_swaps)
    profile = EspressoProfile(argv=ESPRESSO_CMD.split())
    atoms.calc = Espresso(profile = profile,
                          pseudopotentials=pseudopotentials,
                          **flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
