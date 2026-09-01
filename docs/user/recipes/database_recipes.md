# Database-Compatible Recipes

Some materials databases are built from large collections of DFT calculations run with a specific, fixed workflow. To let users generate new entries that are consistent with one of these databases, quacc ships a set of recipes that reproduce those workflows.

## Materials Project

The [Materials Project](https://materialsproject.org/) is a database of inorganic materials. It has used two generations of VASP workflows:

- **Legacy PBE(+U)**: the original PBE(+U) workflow used for the majority of the database.
- **PBEsol/r2SCAN (MP24)**: the newer workflow, which pre-relaxes with PBEsol before a production r2SCAN relaxation and static calculation.

| Workflow | Relax Flow | Static Job | Extra Dependencies |
|---|---|---|---|
| Legacy PBE(+U) | [quacc.recipes.vasp.mp_legacy.mp_gga_relax_flow](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/mp_legacy.html#quacc.recipes.vasp.mp_legacy.mp_gga_relax_flow) | [quacc.recipes.vasp.mp_legacy.mp_gga_static_job](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/mp_legacy.html#quacc.recipes.vasp.mp_legacy.mp_gga_static_job) | `quacc[mp]` |
| PBEsol/r2SCAN (MP24) | [quacc.recipes.vasp.mp24.mp_metagga_relax_flow](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/mp24.html#quacc.recipes.vasp.mp24.mp_metagga_relax_flow) | [quacc.recipes.vasp.mp24.mp_metagga_static_job](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/mp24.html#quacc.recipes.vasp.mp24.mp_metagga_static_job) | `quacc[mp]` |

## QMOF Database

The [QMOF](https://doi.org/10.1016/j.matt.2021.02.015) Database is a collection of metal-organic frameworks (MOFs). [quacc.recipes.vasp.qmof.qmof_relax_job](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/qmof.html#quacc.recipes.vasp.qmof.qmof_relax_job) reproduces the QMOF Database's multi-step relaxation workflow (an optional pre-relaxation, loose position/volume relaxations, a production-quality double relaxation, and a final static calculation) in a single job.

| Database | Relax Flow | Extra Dependencies |
|---|---|---|
| [QMOF](https://doi.org/10.1016/j.matt.2021.02.015) | [quacc.recipes.vasp.qmof.qmof_relax_job](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/qmof.html#quacc.recipes.vasp.qmof.qmof_relax_job) | none |

A standalone QMOF-compatible static calculation can be run through [quacc.recipes.vasp.core.static_job](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/vasp/core.html#quacc.recipes.vasp.core.static_job) directly, with `preset="`[`QMOFSet`](https://github.com/Quantum-Accelerators/quacc/blob/main/src/quacc/calculators/vasp/presets/QMOFSet.yaml)`"`.

## A Representative Example

```python
from ase.build import bulk
from quacc.recipes.vasp.mp24 import mp_metagga_relax_flow

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run the full MP24 (PBEsol/r2SCAN) relaxation workflow
result = mp_metagga_relax_flow(atoms)
print(result)
```
