
## Example 1

```python
from ase.build import bulk
from quacc.recipes.emt.core import static_job

atoms = bulk("Cu")
result = static_job(atoms)
print(result)
```

## Example 2

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

atoms = bulk("Cu")
result1 = relax_job(atoms)
result2 = static_job(result1["atoms"])
print(result2)
```

## Example 3

```python
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

atoms1 = bulk("Cu")
atoms2 = molecule("N2")
result1 = relax_job(atoms1)
result2 = relax_job(atoms2)
print([result1, result2])
```

## Example 4

```python
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_slabs_from_bulk,make_adsorbate_structures

cu = bulk("Cu")
n2 = molecule("N2")
relax_cu = relax_job(cu)
relax_n2 = relax_job(n2)

slabs = make_slabs_from_bulk(relax_cu["atoms"])
results = []
for slab in slabs:
    slab_ads_systems = make_adsorbate_structures(slab,relax_n2["atoms"])
    for slab_ads_system in slab_ads_systems:
        result = static_job(slab_ads_system)
        results.append(result)
print(results)
```
