# MLIP-Compatible Recipes

Several popular datasets used to train machine-learned interatomic potentials (MLIPs) were built from large collections of DFT calculations run with a specific, fixed set of computational settings. To let users generate additional data that is consistent with one of these datasets (e.g. to extend a dataset or build an internally-consistent benchmark without introducing methodological drift), quacc ships a set of recipes that reproduce those settings.

!!! Note

    While these recipes aim to reproduce the settings of each dataset, you should still manually confirm that the resulting calculator settings are consistent with the source dataset for good measure.

## Supported Datasets

| Dataset | Recipe | Extra Dependencies |
|---|---|---|
| [OC20](https://doi.org/10.1021/acscatal.0c04525) | [quacc.recipes.vasp.fairchem.oc20_static_job][] | `quacc[fairchem]` |
| [OC22](https://doi.org/10.1021/acscatal.2c05426) | [quacc.recipes.vasp.fairchem.oc22_static_job][] | none |
| [OC25](https://arxiv.org/abs/2509.17862v1) | [quacc.recipes.vasp.fairchem.oc25_static_job][] | `quacc[fairchem]` |
| [OMat24](https://www.nature.com/articles/s43588-026-00996-w) | [quacc.recipes.vasp.fairchem.omat_static_job][] | `quacc[fairchem]` |
| [OMC25](https://www.nature.com/articles/s41597-026-06628-2) | [quacc.recipes.vasp.fairchem.omc_static_job][] | `quacc[atomate2]` |
| [OMol25](https://arxiv.org/abs/2505.08762) | [quacc.recipes.orca.fairchem.omol_static_job][] | `quacc[fairchem]` |
| [ODAC23](https://doi.org/10.1021/acscentsci.3c01629) | [quacc.recipes.vasp.fairchem.odac23_static_job][] | none |
| [ODAC25](https://arxiv.org/abs/2508.03162) | [quacc.recipes.vasp.fairchem.odac25_static_job][] | none |
| [AQCat25](https://www.nature.com/articles/s41524-026-02099-6) | [quacc.recipes.vasp.aqcat.aqcat25_static_job][] | none |
| [MPtrj](https://www.nature.com/articles/s42256-023-00716-3) / [WBM](https://www.nature.com/articles/s41524-020-00481-6) / [sAlex](https://matbench-discovery.materialsproject.org/data/salex) / [MatterSim](https://arxiv.org/abs/2405.04967) | [quacc.recipes.vasp.mp_legacy.mp_relax_set_job][] | none |
| [MatPES](https://arxiv.org/abs/2503.04070) | [quacc.recipes.vasp.matpes.matpes_static_job][] | `quacc[atomate2]` |
| [MP-ALOE](https://www.nature.com/articles/s41524-025-01834-9) | [quacc.recipes.vasp.mp_aloe.mp_aloe_static_job][] | none |
| [MAD-1.5 / MAD-1.6](https://arxiv.org/abs/2603.02089) | [quacc.recipes.aims.mad.static_job][] | none |

### MAD-1.5 and MAD-1.6

The MAD recipe uses FHI-aims to compute r2SCAN reference energies, forces, and periodic stresses. The published `control.in` files in the [MAD-1.5 archive](https://archive.materialscloud.org/records/18tke-tt476) and [MAD-1.6 archive](https://archive.materialscloud.org/records/ffd5a-1hc37) are identical, so both datasets share this recipe. Structure generation and dataset filtering are outside its scope.

Configure `AIMS_BIN`, `AIMS_PARALLEL_CMD`, and `AIMS_SPECIES_DEFAULTS` as usual. `AIMS_SPECIES_DEFAULTS` should point to the FHI-aims `defaults_2020` directory containing `tight/`. The paper used FHI-aims version 250806; a build with LibXC support for r2SCAN is required. The recipe automatically removes the confined 5d function for Pr-Yb and the confined 6d function for Pu-No using temporary copies of the species defaults, preserving the installed files.

For the `mc2d` and `mc3d_surface` subsets, set `atoms.info["subset"]` to the corresponding name to enable the published dipole correction, or pass `use_dipole_correction=True` explicitly. Slabs and 2D materials need a complete cell containing vacuum. Molecules and isolated atoms omit the k-point density and stress defaults, even if their `Atoms` object has a cell. See [quacc.recipes.aims.mad.static_job][] for the full recipe options.

## A Representative Example

```python
from ase.build import bulk
from quacc.recipes.vasp.fairchem import omat_static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a static calculation with OMat24-compatible VASP settings
result = omat_static_job(atoms)
print(result)
```

As with all quacc recipes, you can override any of the default calculator settings by passing additional keyword arguments:

```python
result = omat_static_job(atoms, encut=600)
```
