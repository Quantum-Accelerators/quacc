# FHI-aims MAD Recipes

The [quacc.recipes.aims.mad.static_job][] recipe reproduces the shared r2SCAN single-point protocol used for the MAD-1.5 and MAD-1.6 datasets. The published `control.in` files in the [MAD-1.5 archive](https://archive.materialscloud.org/records/18tke-tt476) and [MAD-1.6 archive](https://archive.materialscloud.org/records/ffd5a-1hc37) are identical. The recipe computes energies, forces, and stresses for periodic systems. Structure generation and dataset filtering are outside its scope.

Configure `AIMS_BIN`, `AIMS_PARALLEL_CMD`, and `AIMS_SPECIES_DEFAULTS` as usual. `AIMS_SPECIES_DEFAULTS` should point to the FHI-aims `defaults_2020` directory containing `tight/`. The [paper](https://arxiv.org/abs/2603.02089) used FHI-aims version 250806. A build with LibXC support for r2SCAN is required.

```python
from ase.io import read
from quacc.recipes.aims.mad import static_job

atoms = read("geometry.in", format="aims")
result = static_job(atoms)
```

The recipe follows the archived adaptive SCF protocol (`adjust_scf always 2`), including both occupation-width settings of 0.05 eV, non-spin-polarized neutral calculations, scalar atomic ZORA, and the published convergence thresholds. It uses FHI-aims' native `k_grid_density 8` for periodic structures. Molecules and isolated atoms omit the k-point density and stress flags, even if their `Atoms` object has a cell.

For the `mc2d` and `mc3d_surface` subsets, the published input requires a dipole correction. The recipe enables this when `atoms.info["subset"]` has either name. When working with structures without that metadata, enable it explicitly:

```python
result = static_job(slab_atoms, use_dipole_correction=True)
```

Slabs and 2D materials need a complete cell containing vacuum. Geometry, vacuum size, and surface orientation must match the structure being reproduced.

The recipe applies the published species-default modification automatically: it removes the confined 5d function for Pr-Yb and the confined 6d function for Pu-No. Ce and the remaining elements retain their original basis. Modified files are temporary copies; installed FHI-aims species files are preserved. To select another species directory, pass `species_dir` as a path relative to `AIMS_SPECIES_DEFAULTS` or an absolute path.

Calculator overrides, `quacc.Remove`, `copy_files`, and `additional_fields` work as in the core FHI-aims recipes. Overrides can change the level of theory and therefore the agreement with the published dataset.
