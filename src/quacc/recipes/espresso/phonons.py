"""Phonon recipes for espresso."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from quacc import job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import base_fn

if TYPE_CHECKING:
    from quacc.schemas._aliases.ase import RunSchema


@job
def phonon_job(
    prev_dir: str | Path, parallel_info: dict[str] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a basic ph.x calculation.

    Parameters
    ----------
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    **calc_kwargs
        calc_kwargs dictionary possibly containing the following keys:

        - input_data: dict
        - qpts: list[list[float]] | list[tuple[float]] | list[float]
        - nat_todo: list[int]

        See the docstring of ase.io.espresso.write_espresso_ph for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        "input_data": {
            "inputph": {
                "tr2_ph": 1e-12,
                "alpha_mix(1)": 0.1,
                "nmix_ph": 12,
                "verbosity": "high",
            }
        },
        "qpts": (0, 0, 0),
    }

    return base_fn(
        template=EspressoTemplate("ph"),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "ph.x Phonon"},
        copy_files=prev_dir,
    )
