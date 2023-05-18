# Quick Start

Want to get up and running with Quacc as fast possible? Here we go!

## Installation

1. Run `pip install git+https://github.com/arosen93/quacc.git`.
2. Run `covalent start` and open the URL.
3. Run a sample workflow and check the webpage for results.

## Demo Workflow 1: A Simple One

Run the following

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job

# Define the workflow
workflow = ct.lattice(relax_job)

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Dispatch the workflow to the Covalent server
# with the bulk Cu Atoms object as the input
dispatch_id = ct.dispatch(workflow)(atoms)

# Fetch the result from the server
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

## Demo Workflow 2: A More Complex One

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.slabs import BulkToSlabsFlow, relax_job

@ct.lattice(executor="local")
def workflow(atoms):
    relaxed_bulk = relax_job(atoms)
    relaxed_slabs = BulkToSlabsFlow().run(relaxed_bulk.output["atoms"])
    return relaxed_slabs

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
```

## What Next?

Read through the docuemntation to learn more about Quacc and how to use it!