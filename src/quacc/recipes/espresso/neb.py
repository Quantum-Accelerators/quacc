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
from quacc.calculators.espresso.utils import clean_params
from quacc.calculators.espresso.calculator import Neb
from quacc.runners.prep import calc_setup
from ase import Atoms

ESPRESSO_NEB_CMD = f"{SETTINGS.ESPRESSO_NEB_CMD}"
ESPRESSO_PP_PATH = f"{SETTINGS.ESPRESSO_PP_PATH}"
ESPRESSO_PRESET_PATH = f"{SETTINGS.ESPRESSO_PRESET_PATH}"

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def neb_job(
    atoms_list: list[Atoms],
    neb_input_data: dict[str, Any],
    pw_input_data: dict[int, Any],
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
    # List must be size >= 2
    # Preset are assumed to be for all atoms in list

    neb_calc = Neb(neb_input_data,
                   pw_input_data)
    # For some reason this function needs an Atoms object
    # we do that for now
    _, tmpdir, job_results_dir = calc_setup(Atoms(), copy_files=copy_files)

    neb_calc.execute(preset=preset, **kwargs)
    results = neb_calc.read_results(**kwargs)

    calc_cleanup(tmpdir, job_results_dir)
    return _base_job(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=kwargs,
        additional_fields={"name": "pw.x static"},
        copy_files=copy_files,
    )