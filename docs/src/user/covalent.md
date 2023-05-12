# Using Covalent

Here, we will show how to use [Covalent](https://github.com/AgnostiqHQ/covalent) to construct, dispatch, and monitor workflows in QuAcc. This section is intentionally designed to mimic the [Quick Start: Basics](quick-start.md) section so that you can see how it is trivial to go from "normal Python" to a distributed workflow.

## Pre-Requisites

If you haven't already done so, make sure you completed the ["Covalent Setup"](covalent.md) section of the documentation. If you haven't done so already, run `covalent start` prior to starting these examples.

Additionally, you should read the Covalent [First Experiment](https://covalent.readthedocs.io/en/latest/getting_started/first_experiment/index.html) guide to get a sense of how Covalent works. Namely, you should understand the [Covalent Basics](https://covalent.readthedocs.io/en/latest/concepts/basics.html) of the `Electron` and `Lattice` objects, which describe individual compute tasks and workflows, respectively.

## Example 2: A Simple Calculation with GFN2-xTB

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

atoms = bulk("Cu")

@ct.lattice
def workflow(atoms):
    job1 = ct.electron(static_job)
    job2 = ct.electron(relax_job)

    result1 = job1(atoms)
    result2 = job2(result1["atoms"])
    return result2
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

### Example 3: Running a Simple Parallel Workflow

```python
import covalent as ct
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

bulk = bulk("Cu")
molecule = molecule("N2")
@ct.lattice
def workflow(atoms1,atoms2):
    job1 = ct.electron(relax_job)
    job2 = ct.electron(relax_job)

    result1 = job1(atoms1)
    result2 = job2(atoms2)
    return [result1,result2]
dispatch_id = ct.dispatch(workflow)(bulk,molecule)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

### Example 4: Running Workflows with Complex Connectivity

```python
import covalent as ct
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_slabs_from_bulk,make_adsorbate_structures

@ct.lattice
def workflow(bulk,molecule):
    relax_cu = ct.electron(relax_job)(bulk)
    relax_n2 = ct.electron(relax_job)(molecule)

    slabs = ct.electron(make_slabs_from_bulk)(relax_cu["atoms"])
    results = []
    for slab in slabs:
        slab_ads_systems =ct.electron(make_adsorbate_structures)(slab,relax_n2["atoms"])
        for slab_ads_system in slab_ads_systems:
            result = ct.electron(static_job)(slab_ads_system)
            results.append(result)
    return results

cu = bulk("Cu")
n2 = molecule("N2")
dispatch_id = ct.dispatch(workflow)(cu,n2)
results = ct.get_result(dispatch_id, wait=True)
print(results)
```

## Shorthand for Single-Job Lattices

Let's return to the simple example of using EMT to run a static calculation on a bulk structure of copper (Example 1 from the [Quick Start: Basics](quick-start.md) guide). This can be achieved using Covalent as follows:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define an Electron (job) that will perform the static calculation
@ct.electron
def job(atoms):
    result = static_job(atoms)
    return result

# Define a Lattice (workflow) that will consist of the individual job steps
@ct.lattice
def workflow(atoms):
    return job(atoms)

dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

Now, this is admittedly quite verbose. Because this is a very simple one-job workflow, it can be made much more compact as follows:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define a Lattice (workflow) and Electron (job)
workflow = ct.lattice(ct.electron(static_job))
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```
