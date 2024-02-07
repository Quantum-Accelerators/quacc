from ase.units import fs

@job
def microcanonical_flow(
    atoms: Atoms,
    production_job: dict[Any, dict],
    time_step: float = fs,
    initial_temperature: float | None = None,
    equilibration_jobs: dict[Any, dict] | None = None,
) -> MdSchema:
    """
    Perform a microcanonical flow calculation with optional equilibration runs.

    Parameters
    ----------
    atoms
        Atoms object with calculator attached.
    production_job
        The Class to use to run the production job. This should be a dictionary
        where the key is the class and the values are the kwargs to pass to the
        class.
    initial_temperature
        Initial temperature, if None atoms will not have any distribution. In NVE
        it should be set to twice the desired temperature, due to the equipartition theorem.
    equilibration_jobs
        Equilibration jobs to run before the production job. This should be a dictionary
        of the form {job_class: job_kwargs}. The job_kwargs should be a dictionary

    Returns
    -------
    MdSchema
    """

    return None
