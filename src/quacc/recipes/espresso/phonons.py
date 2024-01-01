"""
This module, 'phonons.py', contains recipes for performing phonon calculations using the ph.x binary from Quantum ESPRESSO via the quacc library. The recipes provided in this module are jobs and flows that can be used to perform phonon calculations in different fashion. Please refer to the individual function docstrings for more detailed information on their usage and parameters.

If you don't know how to proceed, please consult [] to learn how to use the
quacc espresso calculator from basics to advanced.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from quacc import Job, flow, job, strip_decorator, subflow
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.calculators.espresso.utils import parse_ph_patterns
from quacc.recipes.espresso._base import base_fn
from quacc.recipes.espresso.core import static_job, relax_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def phonon_job(
    prev_dir: str | Path,
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic ph.x calculation. It should allow you to
    use all the features of the ph.x binary (https://www.quantum-espresso.org/Doc/INPUT_PH.html)

    This job requires the results of a previous pw.x calculation, you might
    want to create your own flow to run both jobs in sequence. If you don't know how to do it please consult [].

    Parameters
    ----------
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
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
        template=EspressoTemplate("ph", test_run=test_run),
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
    This functions is a subflow used in [quacc.recipes.espresso.phonons.grid_phonon_flow][]. Feel free to call it directly in your own flow, for
    example if you already have a pw.x calculation and you simply want to run a grid parallelized ph.x calculation.

    If you don't know how to create your own flow, please consult [].

    Parameters
    ----------
        ph_job:
            The phonon job to be executed.

        pw_job_results_dir
            The directory containing the results of the plane-wave job.

        nblocks
            The number of blocks for grouping representations. Defaults to 1.

    Returns
    -------
        list[RunSchema]: A list of results from each phonon job.
    """
    # This will be forced to run locally, do we want that? FYI, test_run is calling the binaries
    ph_test_job_results = strip_decorator(ph_job)(pw_job_results_dir, test_run=True)
    input_data = ph_test_job_results["parameters"]["input_data"]
    prefix = input_data["inputph"].get("prefix", "pwscf")
    ph_patterns = parse_ph_patterns(ph_test_job_results["dir_name"], prefix)

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
    This function performs grid parallelization of a ph.x calculation. Each representation of each q-point is calculated in a separate job, allowing for distributed computation across different machines and times.

    The grid parallelization is a technique to make phonon calculation embarrassingly parallel. This function should return similar results to [quacc.recipes.espresso.phonons.phonon_job][]. If you don't know about
    grid parallelization please consult the Quantum Espresso user manual and
    exemples.

    This approach requires the data of the pw.x calculation to be copied to each job, leading to a total data size on the disk of n*m times the size of the pw.x calculation, where:
    - n is the number of q-points
    - m is the number of representations

    In addition to the data produced by each ph.x calculation. This can result in large data sizes for systems with many atoms.

    To mitigate this, an optional "nblocks" argument can be provided. This groups multiple representations together in a single job, reducing the data size by a factor of nblocks, but also reducing the level of parallelization. In the case of nblocks = 0, each job will contain all the representations for each q-point.

    Consists of following jobs that can be modified:

    1. pw.x calculation ("pw_job")

    2. ph.x calculation test_run ("ph_job")

    3. (n * m) / nblocks ph.x calculations ("ph_job")

    4. ph.x calculation to gather data and diagonalize each dynamical matrix ("recover_ph_job")

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

    calc_defaults = {
        "relax_job": {
            "input_data": {
                "control": {"forc_conv_thr": 5.0e-5},
                "electrons": {"conv_thr": 1e-12},
            }
        }
    }

    job_params = recursive_dict_merge(calc_defaults, job_params)

    pw_job, ph_job, recover_ph_job = customize_funcs(
        ["pw_job", "ph_job", "recover_ph_job"],
        [relax_job, phonon_job, phonon_job],
        decorators=job_decorators,
        parameters=job_params,
    )

    pw_job_results = pw_job(atoms)

    grid_results = _phonon_subflow(ph_job, pw_job_results["dir_name"], nblocks=nblocks)

    copy_back = [result["dir_name"] for result in grid_results]

    input_data = grid_results[-1]["parameters"]["input_data"]
    for k in ["start_q", "last_q", "start_irr", "last_irr"]:
        input_data["inputph"].pop(k)
    input_data["inputph"]["recover"] = True

    return recover_ph_job(copy_back, input_data=input_data)
