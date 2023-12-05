"""Core recipes for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms

from quacc import SETTINGS, job
from quacc.calculators.espresso.espresso import EspressoProfile, EspressoTemplate
from quacc.recipes.espresso._base import base_fn

if TYPE_CHECKING:
    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = None,
    copy_files: list[str] | None = None,
    parallel_info: dict[str] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **calc_kwargs
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        "input_data": {
            "control": {"calculation": "scf", "restart_mode": "from_scratch"},
            "system": {"ecutwfc": 60, "ecutrho": 240},
            "electrons": {"conv_thr": 1e-8, "mixing_mode": "plain", "mixing_beta": 0.7},
        }
    }

    template = EspressoTemplate("pw")

    return base_fn(
        atoms,
        preset=preset,
        template=template,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x static"},
        copy_files=copy_files,
    )


@job
def ph_job(
    preset: str | None = None,
    copy_files: list[str] | None = None,
    parallel_info: dict[str] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **calc_kwargs
        Custom kwargs for the espresso calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.espresso.Espresso` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    # Default is phonon at gamma, with tight convergence
    # and a lower alphamix(1) as it is very very often
    # recommended in the QE mailing list...
    calc_defaults = {
        "input_data": {
            "inputph": {
                "tr2_ph": 1e-16,
                "alpha_mix(1)": 0.1,
                "nmix_ph": 12,
                "verbosity": "high",
            }
        },
        "qpts": [(0, 0, 0)],
    }

    template = EspressoTemplate("ph")
    profile = EspressoProfile(
        binary=str(SETTINGS.ESPRESSO_PH_CMD),
        parallel_info=parallel_info,
        pseudo_path=None,
    )

    return base_fn(
        Atoms(),
        preset=preset,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "ph.x static"},
        copy_files=copy_files,
    )
