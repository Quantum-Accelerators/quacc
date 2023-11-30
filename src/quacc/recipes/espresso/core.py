"""Core recipes for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.espresso.espresso import Espresso, EspressoTemplate

from quacc import job
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = None,
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
            }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    
    template = EspressoTemplate('pw')

    return _base_job(
        atoms,
        template,
        preset=preset,
        calc_swaps=kwargs,
        additional_fields={"name": "pw.x static"},
        copy_files=copy_files,
    )

@job
def ph_job(
    atoms: Atoms,
    preset: str | None = None,
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
            }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    
    template = EspressoTemplate('ph')

    return _base_job(
        atoms,
        template,
        preset=preset,
        calc_swaps=kwargs,
        additional_fields={"name": "ph.x static"},
        copy_files=copy_files,
    )



def _base_job(
    atoms: Atoms,
    template: EspressoTemplate,
    preset: str | None = None,
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
    # Not happy with this mess currently

    # This is needed because there are two ways
    # of defining the input_data dict, either as a
    # nested dict with each namelist a parent dict or as a flat dict
    # if flat dict, then construct_namelist will convert
    # it to a nested dict, we have to choose one or the other
    # and stick with it, which will be the nested dict here.

    atoms.calc = Espresso(input_atoms = atoms,
                          template = template,
                          preset = preset,
                          **calc_swaps)
    
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
