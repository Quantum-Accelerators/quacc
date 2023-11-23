"""Core recipes for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io.espresso import construct_namelist

from quacc import SETTINGS, job
from quacc.calculators.espresso.io import parse_pp_and_cutoff
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.files import load_yaml_calc

ESPRESSO_CMD = f"{SETTINGS.ESPRESSO_CMD}"
ESPRESSO_PP_PATH = f"{SETTINGS.ESPRESSO_PP_PATH}"
ESPRESSO_PRESET_PATH = f"{SETTINGS.ESPRESSO_PRESET_PATH}"

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
    
    return _base_job(
        atoms,
        preset=preset,
        calc_swaps=kwargs,
        additional_fields={"name": "pw.x static"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms,
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
    input_data = calc_swaps['input_data']
    input_data = construct_namelist(input_data)
    input_data['control']['pseudo_dir'] = ESPRESSO_PP_PATH

    if preset:
        config = load_yaml_calc(ESPRESSO_PRESET_PATH / f"{preset}")
        preset_pp = parse_pp_and_cutoff(config, atoms)
        input_data = input_data | preset_pp['input_data']
        calc_swaps = calc_swaps | preset_pp['pseudopotentials']

    calc_swaps["input_data"] = input_data

    profile = EspressoProfile(argv=ESPRESSO_CMD.split())
    atoms.calc = Espresso(profile = profile,
                          **calc_swaps)
    
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
