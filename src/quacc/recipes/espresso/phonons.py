"""Phonon recipes for espresso."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from quacc import Job, flow, job, strip_decorator, subflow
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.calculators.espresso.utils import parse_ph_patterns
from quacc.recipes.espresso._base import base_fn
from quacc.recipes.espresso.core import static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def phonon_job(
    prev_dir: str | Path,
    preset: str | None = "basic",
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic ph.x calculation.

    Parameters
    ----------
    prev_dir
        Previous directory where pw.x was run. This is used to copy the entire tree structure
        of that directory to the working directory of this calculation.
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    **calc_kwargs
        calc_kwargs dictionary possibly containing the following keys:

        - input_data: dict
        - qpts: list[list[float]] | list[tuple[float]] | list[float]
        - nat_todo: list[int]

        See the docstring of [quacc.calculators.espresso.io.write_espresso_ph][] for more information.

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

    template = EspressoTemplate("ph", test_run=test_run)

    return base_fn(
        preset=preset,
        template=template,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "ph.x Phonon"},
        copy_files=prev_dir,
    )


@subflow
def _phonon_subflow(
    ph_job: Job, pw_job_results_dir: str | Path, nblocks: int = 1
) -> list[RunSchema]:
    """
    TODO.
    """
    # We loop over the q-points and representations and run a ph job for each
    # of them

    # Run a test phonon job
    ph_test_job_results = strip_decorator(ph_job)(pw_job_results_dir, test_run=True)
    input_data = ph_test_job_results["parameters"]["input_data"]
    ph_patterns = parse_ph_patterns(ph_test_job_results["dir_name"])

    grid_results = []
    for pattern in ph_patterns:
        input_data["inputph"]["start_q"] = pattern
        input_data["inputph"]["last_q"] = pattern
        n_repr = ph_patterns[pattern]
        this_block = nblocks if nblocks > 0 else n_repr
        repr_to_do = np.array_split(
            np.arange(1, n_repr + 1), np.ceil(n_repr / this_block)
        )
        for representation in repr_to_do:
            input_data["inputph"]["start_irr"] = representation[0]
            input_data["inputph"]["last_irr"] = representation[-1]
            ph_job_results = ph_job(pw_job_results_dir, input_data=input_data)
            grid_results.append(ph_job_results)
    return grid_results


@flow
def grid_phonon_flow(
    atoms: Atoms,
    nblocks: int = 1,
    job_decorators: dict[str, Callable | None] | None = None,
    job_params: dict[str, Any] | None = None,
) -> RunSchema:
    """
    Function to carry out a grid parallelization of a ph.x calculation. Each representation of each
    q-point is calculated in a separate job. This allow the calculations to be run on different machines
    at possibly different times. The only drawback is that the data of the prior pw.x calculation has to
    be copied to each of the jobs. Thus the total amount of data will be n*m times the size of the pw.x
    calculation where n is the number of q-points and m is the number of representations. Plus the data
    produced by each ph.x calculations. This can becomes very large if your system has a lot of atoms.
    To provide a way to control the amount of data produced, the user can provide an additional argument
    "nblocks" which will group multiple representations together in a single job. This will reduce the
    amount of data produced by a factor of nblocks.

    Consists of following jobs that can be modified:

    1. pw.x calculation ("pw_job")

    2. ph.x calculation ("ph_job")

    3. ph.x calculation to diagonalize the dynamical matrix ("recover_ph_job")

    Parameters
    ----------
    atoms
        Atoms object
    nblocks
        The number of representations to group together in a single job. This will reduce the amount
        of data produced by a factor of nblocks. If nblocks = 0 each job will contain all the representations
        for a single q-point.
    job_decorators
        Custom decorators to apply to each Job in the Flow.
        Refer to [quacc.wflow_tools.customizers.customize_funcs][] for details.
    job_params
        Custom parameters to pass to each Job in the Flow.
        Refer to [quacc.wflow_tools.customizers.customize_funcs][] for details.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    pw_job, ph_job, recover_ph_job = customize_funcs(
        ["pw_job", "ph_job", "recover_ph_job"],
        [static_job, phonon_job, phonon_job],
        decorators=job_decorators,
        parameters=job_params,
    )

    # We run a first pw job (single point or relax) depending on the input
    # ASR: Where is the relax job??
    pw_job_results = pw_job(atoms)

    # Run the phonon subflow
    grid_results = _phonon_subflow(ph_job, pw_job_results["dir_name"], nblocks=nblocks)

    # We have to copy back the files from all of the jobs
    copy_back = [result["dir_name"] for result in grid_results]

    # We run a 'recover' ph job to diagonalize the dynamical matrix
    # for each q-point.
    input_data = grid_results[-1]["parameters"]["input_data"]
    for k in ["start_q", "last_q", "start_irr", "last_irr"]:
        input_data["inputph"].pop(k)
    input_data["inputph"]["recover"] = True

    return recover_ph_job(copy_back, input_data=input_data)
