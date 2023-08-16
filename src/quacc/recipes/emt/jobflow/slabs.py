"""Slab recipes for EMT based on Jobflow"""
from __future__ import annotations

from ase import Atoms
from jobflow import Flow, Job, Response, job

from quacc.recipes.emt.core import relax_job, static_job
from quacc.schemas.atoms import fetch_atoms
from quacc.util.slabs import make_max_slabs_from_bulk


def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    make_slabs_kwargs: dict | None = None,
    slab_relax: Job = job(relax_job),
    slab_static: Job | None = job(static_job),
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
    make_slabs_kwargs
        Additional keyword arguments to pass to `make_max_slabs_from_bulk()`
    slab_relax
        Maker to use for the relaxation of the slab.
    slab_static
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
    atoms = fetch_atoms(atoms)
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    # Generate all the slab
    slabs = make_max_slabs_from_bulk(atoms, **make_slabs_kwargs)

    # Generate the jobs for each slab
    jobs = []
    outputs = []
    for slab in slabs:
        if slab_static is None:
            job1 = slab_relax(slab, **slab_relax_kwargs)
            jobs += [job1]
            outputs.append(job1.output)
        else:
            job1 = slab_relax(slab, **slab_relax_kwargs)
            job2 = slab_static(job1.output, **slab_static_kwargs)
            jobs += [job1, job2]
            outputs.append(job2.output)

    return Response(
        output={"input_bulk": atoms, "generated_slabs": slabs},
        replace=Flow(jobs, output=outputs),
    )
