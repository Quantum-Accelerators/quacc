"""Slab recipes for EMT based on Jobflow"""
from __future__ import annotations

from ase import Atoms
from jobflow import Flow, Response, job

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_max_slabs_from_bulk


def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    slabgen_kwargs: dict | None = None,
    slab_relax_job: job = job(relax_job),
    slab_static_job: job | None = job(static_job),
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> Response:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    slabgen_kwargs
        Additional keyword arguments to pass to `make_max_slabs_from_bulk()`
    slab_relax_job
        Maker to use for the relaxation of the slab.
    slab_static_job
        Maker to use for the static calculation of the slab.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    Response
        A Response containing Flow of relaxation and static jobs for the generated slabs.
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    slabgen_kwargs = slabgen_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    # Generate all the slab
    slabs = make_max_slabs_from_bulk(atoms, **slabgen_kwargs)

    # Generate the jobs for each slab
    jobs = []
    outputs = []
    for slab in slabs:
        if slab_static_job is None:
            job1 = slab_relax_job(slab, **slab_relax_kwargs)
            jobs += [job1]
            outputs.append(job1.output)
        else:
            job1 = slab_relax_job(slab, **slab_relax_kwargs)
            job2 = slab_static_job(job1.output, **slab_static_kwargs)
            jobs += [job1, job2]
            outputs.append(job2.output)

    return Response(
        output={"input_bulk": atoms, "generated_slabs": slabs},
        replace=Flow(jobs, output=outputs),
    )
