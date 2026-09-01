# MLIP-Compatible Recipes

Several popular datasets used to train machine-learned interatomic potentials (MLIPs) were built from large collections of DFT calculations run with a specific, fixed set of computational settings. To let users generate additional data that is consistent with one of these datasets (e.g. to extend a dataset or build an internally-consistent benchmark without introducing methodological drift), quacc ships a set of recipes that reproduce those settings.

## Supported Datasets

| Dataset | Recipe | Extra Dependencies |
|---|---|---|
| [OC20](https://pubs.acs.org/accacs/article/11/10/6059/459430/Open-Catalyst-2020-OC20-Dataset-and-Community) | quacc.recipes.vasp.fairchem.oc20_static_job | `quacc[fairchem]` |
| [OMat24](https://www.nature.com/articles/s43588-026-00996-w) | quacc.recipes.vasp.fairchem.omat_static_job | `quacc[fairchem]` |
| [OMC25](https://www.nature.com/articles/s41597-026-06628-2) | quacc.recipes.vasp.fairchem.omc_static_job | `quacc[mp]` |
| [OMol25](https://arxiv.org/abs/2505.08762) | quacc.recipes.orca.fairchem.omol_static_job | `quacc[fairchem]` |
| [ODAC25](https://arxiv.org/abs/2508.03162) | quacc.recipes.vasp.fairchem.odac_static_job | none |
| [MPtrj](https://www.nature.com/articles/s42256-023-00716-3) / [WBM](https://www.nature.com/articles/s41524-020-00481-6) / [sAlex](https://matbench-discovery.materialsproject.org/data/salex) / [MatterSim](https://arxiv.org/abs/2405.04967) | quacc.recipes.vasp.mp_legacy.mp_relax_set_job | none |
| [MatPES](https://arxiv.org/abs/2503.04070) | quacc.recipes.vasp.matpes.matpes_static_job | `quacc[mp]` |
| [MP-ALOE](https://www.nature.com/articles/s41524-025-01834-9) | quacc.recipes.vasp.mp_aloe.mp_aloe_static_job | none |

## Custodian and the INCAR Copilot

All of the VASP recipes above run with `incar_copilot_mode="critical"` by default, unlike quacc's other VASP recipes. This disables quacc's automatic INCAR-fixing behavior so that the settings intentionally chosen to match a dataset are not silently overridden.

A subset of these recipes — [quacc.recipes.vasp.fairchem.oc20_static_job][], [quacc.recipes.vasp.fairchem.omc_static_job][], and [quacc.recipes.vasp.fairchem.odac_static_job][] — additionally run with `use_custodian=False` by default, since Custodian-based error correction was not used when those datasets were originally generated. The other VASP recipes in the table above leave Custodian enabled (the quacc default).

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
