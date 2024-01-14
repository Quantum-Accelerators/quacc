from pathlib import Path

from quacc.recipes.espresso.core import post_processing_job
from quacc.utils.dicts import recursive_dict_merge


class PostProcessingPatch:
    job_names = ["post_processing"]

    def __init__(self, **calc_kwargs):
        self._results = {job_name: {} for job_name in self.job_names}
        self._calc_kwargs = calc_kwargs

    def step(self, *args, **kwargs):
        super().step(*args, **kwargs)
        # Would this cause issues because everything is ran asynchrously?
        # If the opt is fast compared to the post-processing, file could be
        # overwritten while being copied or something like that?
        # How to redecorate here?
        post_processing_job_params = self._calc_kwargs.get("post_processing_params", {})
        forced_params = {
            "input_data": {
                "plot": {"fileout": f"pseudo_charge_density_{self.nsteps}.cube"}
            }
        }
        self._calc_kwargs["post_processing_params"] = recursive_dict_merge(
            post_processing_job_params, forced_params
        )
        self._results["post_processing"][self.nsteps] = post_processing_job(
            prev_dir=Path().cwd(), **self._calc_kwargs.get("post_processing_params", {})
        )


class SchemaPatch:

    """Base class to patch an ASE optimizer. By modifying methods below,
    you can 'patch' an optimizer to perform custom actions at specific moment of
    the simulation"""

    job_names = []

    def __init__(self, **calc_kwargs):
        self._results = {job_name: {} for job_name in self.job_names}
        self._calc_kwargs = calc_kwargs

    def step(self, *args, **kwargs):
        super().step(*args, **kwargs)

    def run(self, *args, **kwargs):
        super().run(*args, **kwargs)

    def converged(self, *args, **kwargs):
        super().converged(*args, **kwargs)
