# MLIP-Compatible Recipes

Several popular datasets used to train machine-learned interatomic potentials (MLIPs) were built from large collections of DFT calculations run with a specific, fixed set of computational settings. To let users generate additional data that is consistent with one of these datasets (e.g. to extend a dataset or build an internally-consistent benchmark without introducing methodological drift), quacc ships a set of recipes that reproduce those settings.

## Supported Datasets

| Dataset | Recipe |
|---|---|
| OC20 | quacc.recipes.vasp.fairchem.oc20_static_job |
| OMat24 | quacc.recipes.vasp.fairchem.omat_static_job |
| OMC25 | quacc.recipes.vasp.fairchem.omc_static_job |
| OMol25 | quacc.recipes.orca.fairchem.omol_static_job |
| ODAC25 | quacc.recipes.vasp.fairchem.odac_static_job |
| MPtrj | quacc.recipes.vasp.mp_legacy.mp_gga_static_job |
| MatPES | quacc.recipes.vasp.matpes.matpes_static_job |
| MP-ALOE | quacc.recipes.vasp.mp_aloe.mp_aloe_static_job |

A few notes on specific rows:

- **OMol25** is not a VASP recipe. OMol is a molecular (non-periodic) dataset computed with ORCA, so [quacc.recipes.orca.fairchem.omol_static_job][] wraps an ORCA single-point calculation rather than a VASP input set. It requires the `fairchem-data-omol` package (`pip install quacc[fairchem]`).

- **ODAC25**: [quacc.recipes.vasp.fairchem.odac_static_job][]'s default settings were verified against the official `fairchem-data-odac` package and the ODAC25 paper's supplementary Table S2. They match closely, with one known discrepancy (`NELM=60` in quacc vs. `NELM=120` in the paper). Separately, the ODAC25 paper applies a post-hoc statistical correction for k-point convergence to its raw DFT energies — this correction is not a VASP setting and cannot be reproduced by any VASP recipe.

- **MPtrj** does not have a dataset-specific recipe name. [quacc.recipes.vasp.mp_legacy.mp_gga_static_job][] is quacc's existing Materials Project legacy GGA static recipe, which its docstring notes is "also the settings compatible with the MPtrj dataset."

## Custodian and the INCAR Copilot

All of the VASP recipes above run with `incar_copilot_mode="critical"` by default, unlike quacc's other VASP recipes. This disables quacc's automatic INCAR-fixing behavior so that the settings intentionally chosen to match a dataset are not silently overridden.

A subset of these recipes — [quacc.recipes.vasp.fairchem.oc20_static_job][], [quacc.recipes.vasp.fairchem.omc_static_job][], and [quacc.recipes.vasp.fairchem.odac_static_job][] — additionally run with `use_custodian=False` by default, since Custodian-based error correction was not used when those datasets were originally generated. The other VASP recipes in the table above leave Custodian enabled (the quacc default).

## A Representative Example

```python
from ase.build import bulk
from quacc.recipes.vasp.fairchem import oc20_static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a static calculation with OC20-compatible VASP settings
result = oc20_static_job(atoms)
print(result)
```

As with all quacc recipes, you can override any of the default calculator settings by passing additional keyword arguments:

```python
result = oc20_static_job(atoms, encut=600)
```
