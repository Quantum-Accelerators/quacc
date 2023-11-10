# Flow Complexities

## Introduction

There are some important subtleties that appear when writing complex flows. Here, we will consider a simple toy example to demonstrate common pitfalls.

Specifically, let's consider the following example where we:

1. Make many different surface slabs given a bulk material.
2. Run a structure relaxation on each surface slab.

We will call our workflow `bulk_to_slabs_flow`. The convention in quacc is such that the suffix (i.e. `_job` or `_flow`) indicates the type of function it is.

For simplicity, let's assume we already have a `relax_job` and `static_job` defined elsewhere as well as a pre-made function `make_slabs_from_bulk` that will generate the slabs. In other words, let's focus only on the orchestration logic here.

## Correct Attempt

The following is one successful implementation of this toy recipe, which will serve as a useful point of reference.

```python
from quacc import flow
from quacc.atoms.slabs import make_slabs_from_bulk
from quacc.recipes.emt.core import relax_job

@flow
def bulk_to_slabs_flow(atoms):

    @subflow
    def _relax_job_distributed(atoms):
        slabs = make_slabs_from_bulk(atoms)
        return [relax_job(slab) for slab in slabs]

    return _relax_job_distributed(atoms)
```

In the above `bulk_to_slabs_flow`, we have a flow that takes an input `Atoms` object, makes a bunch of slabs (i.e. `slabs: list[Atoms]`), and each individual `slab` is relaxed via the `relax_job`. The output of the parent flow is a list of schemas: one for each slab relaxation.

## Pitfall 1: Lack of Concurency

```python
from quacc import flow
from quacc.atoms.slabs import make_slabs_from_bulk
from quacc.recipes.emt.core import relax_job

@subflow
def bulk_to_slabs_flow(atoms):

    slabs = make_slabs_from_bulk(atoms)

    outputs = []
    for slab in slabs:
        outputs.append(relax_job(slab))

    return outputs
```

Arguably the simplest way one might think to write the flow is shown above. This approach has a few problems, but the one we want to highlight here has to do with
