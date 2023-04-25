# Covalent

## Example 1: Running a Single Function

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import static_job

atoms = bulk("Cu")

@ct.lattice
@ct.electron
def workflow(atoms):
    results = static_job(atoms)
    return results

dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

Note that this is can be made more compact as follows:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import static_job

atoms = bulk("Cu")

job = ct.electron(static_job)
workflow = ct.lattice(job)
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

## Example 2: Running a Simple Serial Workflow

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

## Storing the Results in a Database

Covalent automatically stores all the inputs and outputs in an SQLite database, which you can find at the `"db_path"` when you run `covalent config`. However, if you want to store the results in a different database of your choosing, you can use [`maggma`](https://github.com/materialsproject/maggma) to do so quite easily. This is the generally recommended approach. An example is shown below for storing the results in a MongoDB. For assistance with setting up a MongoDB of your own, please see <>.

```python
import covalent as ct
from maggma.stores import MongoStore

# Connect to the database
database = "my_db"
collection_name = "my_collection"
store = MongoStore(database,collection_nam,host="localhost",port=27017,username="my_username",password="my_password")
store.connect()

# Fetch the results
results_dir = ct.get_config()["dispatcher"]["results_dir"]
docs = []
for dispatcher_id in os.listdir(results_dir):
    result = ct.get_result(dispatch_id).result
    docs.append({"dispatcher_id":dispatcher_id, "result": result})

# Store the results
with store:
    store.update(docs, key="dispatcher_id")

# Close the database connection
store.close()
```
