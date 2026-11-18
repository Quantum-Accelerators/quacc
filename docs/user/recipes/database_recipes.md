# Database-Compatible Recipes

Some materials databases are built from large collections of DFT calculations run with a specific, fixed workflow. To let users generate new entries that are consistent with one of these databases, quacc ships a set of recipes that reproduce those workflows.

## Materials Project

The [Materials Project](https://materialsproject.org/) (MP) is a database of computed properties for inorganic materials. Over the years, MP has computed its data with several generations of VASP settings, and quacc provides a recipe for each of them. To combine your own calculations with MP data (e.g. to compute an energy above the convex hull), your calculations must use the same settings, including the same pseudopotentials, as the MP entries you compare against.

| Workflow | Relax Flow | Static Job | Input Sets (pymatgen) | Pseudopotentials |
|---|---|---|---|---|
| Legacy PBE(+U) | [quacc.recipes.vasp.mp_legacy.mp_gga_relax_flow][] | [quacc.recipes.vasp.mp_legacy.mp_gga_static_job][] | `MPRelaxSet`, `MPStaticSet` | Original PBE (`potpaw_PBE.original`) |
| Legacy r2SCAN | [quacc.recipes.vasp.mp_legacy.mp_metagga_relax_flow][] | [quacc.recipes.vasp.mp_legacy.mp_metagga_static_job][] | `MPScanRelaxSet`, `MPScanStaticSet` | PBE.54 (`potpaw_PBE.54`) |
| MP24 r2SCAN | [quacc.recipes.vasp.mp24.mp_metagga_relax_flow][] | [quacc.recipes.vasp.mp24.mp_metagga_static_job][] | `MP24RelaxSet`, `MP24StaticSet` | PBE.64 (`potpaw_PBE.64`) |

All of these recipes require `pip install quacc[atomate2]`.

- **Legacy PBE(+U)**: the original PBE and PBE+U workflow used for most of the GGA data in MP. The flow runs two PBE(+U) relaxations followed by a static calculation, matching atomate2's `MPGGADoubleRelaxStaticMaker`.

- **Legacy r2SCAN**: the r2SCAN workflow used to compute the r2SCAN data currently in MP's core database. The flow runs a PBEsol pre-relaxation, an r2SCAN relaxation, and an r2SCAN static calculation, matching atomate2's `MPMetaGGADoubleRelaxStaticMaker`.

- **MP24 r2SCAN**: the r2SCAN settings MP adopted after its 2023–2024 benchmarking effort. The flow runs a PBEsol pre-relaxation, an r2SCAN relaxation, and an r2SCAN static calculation, matching atomate2's `MP24DoubleRelaxStaticMaker`.

!!! Tip "Which r2SCAN Workflow Should I Use?"

    To compare against or combine with the r2SCAN entries currently in the Materials Project database (e.g. to compute an energy above the convex hull), use the legacy r2SCAN recipes until there has been a complete transition to the MP24 settings from the Materials Project.

!!! Note "Related Datasets"

    Recipes for the MatPES and MP-ALOE datasets, which are also generated with MP-style settings, are listed on the [MLIP recipes](mlip_recipes.md) page.

## QMOF Database

The [QMOF](https://doi.org/10.1016/j.matt.2021.02.015) Database is a collection of computed properties for metal-organic frameworks (MOFs). [quacc.recipes.vasp.qmof.qmof_relax_job][] reproduces the QMOF Database's multi-step relaxation workflow (an optional pre-relaxation, loose position/volume relaxations, a production-quality double relaxation, and a final static calculation) in a single job.

| Database | Relax Flow | Extra Dependencies |
|---|---|---|
| [QMOF](https://doi.org/10.1016/j.matt.2021.02.015) | [quacc.recipes.vasp.qmof.qmof_relax_job][] | none |

A standalone QMOF-compatible static calculation can be run through [quacc.recipes.vasp.core.static_job][] directly, with `preset="QMOFSet"` (see [QMOFSet.yaml](https://github.com/Quantum-Accelerators/quacc/blob/main/src/quacc/calculators/vasp/presets/QMOFSet.yaml)).

## A Representative Example

```python
from ase.build import bulk
from quacc.recipes.vasp.mp_legacy import mp_metagga_relax_flow

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run the legacy r2SCAN workflow, which is compatible with the
# r2SCAN data currently in the Materials Project database
result = mp_metagga_relax_flow(atoms)

# The final energy comes from the r2SCAN static calculation
print(result["static"]["results"]["energy"])
```

To use the MP24 settings instead, import `mp_metagga_relax_flow` from `quacc.recipes.vasp.mp24`. Both flows return a dictionary with `prerelax`, `relax`, and `static` keys.
