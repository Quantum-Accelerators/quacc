"""Core recipes for espresso."""
from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

import numpy as np

from quacc import Job, flow, job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.calculators.espresso.utils import parse_ph_patterns
from quacc.recipes.espresso._base import base_fn
from quacc.recipes.espresso.core import static_job

if TYPE_CHECKING:
    from quacc.schemas._aliases.ase import RunSchema


@job
def phonon_job(
    preset: str | None = None,
    copy_files: str | list[str] | None = None,
    parallel_info: dict[str] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic ph.x calculation.

    Parameters
    ----------
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
    copy_files
        List of files to copy to the calculation directory. Almost always needed
        for ph.x calculations. This parameter can either be a string or a list of
        strings.

        If a string is provided, it is assumed to be a path to a directory,
        all of the child tree structure of that directory is going to be copied to the
        scratch of this calculation. For phonon_job this is what most users will want to do.

        If a list of strings is provided, each string point to a specific file. In this case
        it is important to note that no directory structure is going to be copied, everything
        is copied at the root of the temporary directory.
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
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
        copy_files=copy_files,
    )


@flow
def grid_phonon_flow(
    custom_pw_job: Job | None = None,
    custom_phonon_job: Job | None = None,
    nblocks: int = 1,
):
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

    Parameters
    ----------
    custom_pw_job
        A custom job to run the pw.x calculation. This job should return a RunSchema dictionary.
        The default job is [quacc.recipes.espresso.core.static_job][].
    custom_phonon_job
        A custom job to run the ph.x calculation. This job should return a RunSchema dictionary.
        The default job is [quacc.recipes.espresso.phonons.phonon_job][].
    nblocks
        The number of representations to group together in a single job. This will reduce the amount
        of data produced by a factor of nblocks. If nblocks = 0 each job will contain all the representations
        for a single q-point.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    # We run a first pw job (single point or relax) depending on the input
    pw_job = static_job if custom_pw_job is None else custom_phonon_job
    ph_job = phonon_job if custom_phonon_job is None else custom_phonon_job

    pw_job_results = pw_job()

    # First ph job to generate the patterns
    ph_job_partial = partial(
        ph_job, test_run=True, copy_files=pw_job_results["dir_name"]
    )
    # We actually run the job
    ph_job_results = ph_job_partial()
    # We parse the patterns
    ph_patterns = parse_ph_patterns(ph_job_results["dir_name"])
    # We need to get the input_data from the ph job that was sent in...
    input_data = ph_job_results["parameters"]["input_data"]  # Better way to do that?

    grid_results = []
    # We loop over the q-points and representations and run a ph job for each
    # of them
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
            ph_job_partial = partial(
                ph_job, input_data=input_data, copy_files=pw_job_results["dir_name"]
            )
            grid_results.append(ph_job_partial())
    # We have to copy back the files from all of the jobs
    copy_back = [result["dir_name"] for result in grid_results]
    # We run a 'recover' ph job to diagonalize the dynamical matrix
    # for each q-point.
    input_data["inputph"].pop("start_q")
    input_data["inputph"].pop("last_q")
    input_data["inputph"].pop("start_irr")
    input_data["inputph"].pop("last_irr")

    input_data["inputph"]["recover"] = True

    ph_job_partial = partial(ph_job, input_data=input_data, copy_files=copy_back)

    return ph_job_partial()
