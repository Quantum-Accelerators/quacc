"""Slab recipes for EMT based on Jobflow"""
from __future__ import annotations

from dataclasses import dataclass

from ase import Atoms
from jobflow import Flow, Maker, Response, job

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_max_slabs_from_bulk


@dataclass
class BulkToSlabsFlow(Maker):
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    name
        Name of the job.
    slab_relax_job
        Maker to use for the relaxation of the slab.
    slab_static_job
        Maker to use for the static calculation of the slab.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    name: str = "EMT BulkToSlabsFlow"
    slab_relax_job: job = job(relax_job)
    slab_static_job: job | None = job(static_job)
    slab_relax_kwargs: dict | None = None
    slab_static_kwargs: dict | None = None

    @job
    def make(self, atoms: Atoms, slabgen_kwargs: dict = None) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            Atoms object
        slabgen_kwargs
            Additional keyword arguments to pass to `make_max_slabs_from_bulk()`

        Returns
        -------
        Response
            A Response containing Flow of relaxation and static jobs for the generated slabs.
        """
        self.slab_relax_kwargs = self.slab_relax_kwargs or {}
        self.slab_static_kwargs = self.slab_static_kwargs or {}
        slabgen_kwargs = slabgen_kwargs or {}

        if "relax_cell" not in self.slab_relax_kwargs:
            self.slab_relax_kwargs["relax_cell"] = False

        # Generate all the slab
        slabs = make_max_slabs_from_bulk(atoms, **slabgen_kwargs)

        # Generate the jobs for each slab
        jobs = []
        outputs = []
        for slab in slabs:
            if self.slab_static_job is None:
                job1 = self.slab_relax_job(slab, **self.slab_relax_kwargs)
                jobs += [job1]
                outputs.append(job1.output)
            else:
                job1 = self.slab_relax_job(slab, **self.slab_relax_kwargs)
                job2 = self.slab_static_job(
                    job1.output["atoms"], **self.slab_static_kwargs
                )
                jobs += [job1, job2]
                outputs.append(job2.output)

        return Response(
            output={"input_bulk": atoms, "generated_slabs": slabs},
            replace=Flow(
                jobs,
                output=outputs,
                name=self.name,
            ),
        )
