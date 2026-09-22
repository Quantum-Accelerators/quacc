# MLIP-Compatible Recipes

Several popular datasets used to train machine-learned interatomic potentials (MLIPs) were built from large collections of DFT calculations run with a specific, fixed set of computational settings. To let users generate additional data that is consistent with one of these datasets (e.g. to extend a dataset or build an internally-consistent benchmark without introducing methodological drift), quacc ships a set of recipes that reproduce those settings.

## Supported Datasets

| Dataset | Recipe | Extra Dependencies |
|---|---|---|
| [OC20](https://pubs.acs.org/accacs/article/11/10/6059/459430/Open-Catalyst-2020-OC20-Dataset-and-Community) | [quacc.recipes.vasp.fairchem.oc20_static_job][] | `quacc[fairchem]` |
| [OC25](https://arxiv.org/abs/2509.17862) | [quacc.recipes.vasp.fairchem.oc25_static_job][] | `quacc[fairchem]` |
| [OMat24](https://www.nature.com/articles/s43588-026-00996-w) | [quacc.recipes.vasp.fairchem.omat_static_job][] | `quacc[fairchem]` |
| [OMC25](https://www.nature.com/articles/s41597-026-06628-2) | [quacc.recipes.vasp.fairchem.omc_static_job][] | `quacc[fairchem]` |
| [OMol25](https://arxiv.org/abs/2505.08762) | [quacc.recipes.orca.fairchem.omol_static_job][] | `quacc[fairchem]` |
| [ODAC25](https://arxiv.org/abs/2508.03162) | [quacc.recipes.vasp.fairchem.odac_static_job][] | none |
| [AQCat25](https://www.nature.com/articles/s41524-026-02099-6) | [quacc.recipes.vasp.aqcat.aqcat25_static_job][] | none |
| [MPtrj](https://www.nature.com/articles/s42256-023-00716-3) / [WBM](https://www.nature.com/articles/s41524-020-00481-6) / [sAlex](https://matbench-discovery.materialsproject.org/data/salex) / [MatterSim](https://arxiv.org/abs/2405.04967) | [quacc.recipes.vasp.mp_legacy.mp_relax_set_job][] | none |
| [MatPES](https://arxiv.org/abs/2503.04070) | [quacc.recipes.vasp.matpes.matpes_static_job][] | `quacc[atomate2]` |
| [MP-ALOE](https://www.nature.com/articles/s41524-025-01834-9) | [quacc.recipes.vasp.mp_aloe.mp_aloe_static_job][] | none |

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

For an OC25 interface, use `oc25_static_job(atoms)`. This recipe uses RPBE with
D3 zero damping, the OC25 surface k-point mesh and dipole correction, and ASE's
recommended VASP 6.4 potentials. Configure `VASP_PP_PATH` with a
`potpaw_PBE.64` directory as described in the [VASP installation guide](../../install/codes.md).
Its default `ediff=1e-6` matches the OC25 validation/test single points; passing
`ediff=1e-4` uses the training data's electronic threshold while keeping a static
calculation. Automatic INCAR adjustments and Custodian corrections are disabled
by default to preserve these settings.

## Using OC25-Trained Potentials

The VASP recipes above generate DFT reference data. To predict energies and forces with an OC25-trained potential, use the existing [MLIP static and relaxation recipes][quacc.recipes.mlip.core] and [MD recipe][quacc.recipes.mlip.md.md_job] with `library="fairchem"`.

Install `quacc[fairchem]` and obtain access to the [OC25 eSEN checkpoints](https://huggingface.co/facebook/OC25) or [UMA checkpoints](https://huggingface.co/facebook/UMA). Configure Hugging Face authentication on the machine that will execute the jobs. The models run on CUDA or CPU; CPU execution can be slow.

The conservative OC25 eSEN checkpoint has a single task, selected automatically during local inference. UMA requires `task_name="oc25"` and a checkpoint that includes that task, such as `uma-s-1p2`.

```python
from ase.io import read
from quacc.recipes.mlip.core import relax_job, static_job
from quacc.recipes.mlip.md import md_job

# An ASE-readable structure exported from OC25, including its constraints.
atoms = read("oc25_structure.traj")

model = {
    "library": "fairchem",
    "name_or_path": "esen-sm-conserving-all-oc25",
    "device": "cuda",  # Use "cpu" when a CUDA device is unavailable.
}
# To use UMA instead, replace model with:
# model = {
#     "library": "fairchem",
#     "name_or_path": "uma-s-1p2",
#     "task_name": "oc25",
#     "device": "cuda",
# }

provenance = {"oc25_model": model["name_or_path"]}
static_result = static_job(atoms, additional_fields=provenance, **model)
relax_result = relax_job(
    atoms,
    relax_cell=False,
    opt_params={"fmax": 0.05, "max_steps": 200},
    additional_fields=provenance,
    **model,
)
md_result = md_job(
    relax_result["atoms"],
    dynamics="nvt_langevin",
    steps=100,
    timestep_fs=0.5,
    temperature_K=300,
    additional_fields=provenance,
    **model,
)
```

Retain the source cell, periodic boundary conditions, atom order, tags, and fixed-atom constraints when loading structures. Keep reference DFT energies and forces separately from predictions, and record the source structure identifier and checkpoint version for reproducibility. The recipes copy their input atoms and preserve ASE constraints. Their force results are raw calculator forces; ASE applies constraints during relaxation and dynamics.

These examples keep the cell fixed. OC25 models were not trained on stress labels, so cell relaxation and NPT require separate validation. UMA's OC25 task expects zero charge and spin model inputs; it does not predict work functions. See the [FAIR-Chem model documentation](https://github.com/facebookresearch/fairchem/blob/main/docs/core/uma.md) for the training domain and limitations. FAIR-Chem's `free_energy` is an alias for predicted energy, not a separately predicted electronic free energy.

For local inference, quacc caches calculators by their arguments within each worker process. A relaxation or MD trajectory reuses its calculator at every step. This does not share a model across separate worker processes; account for one model copy per worker when assigning GPU resources.

### Optional Integration Tests

The OC25 integration tests compare direct FAIR-Chem and quacc predictions on the same structure, then exercise constrained relaxation and short fixed-cell MD for conservative eSEN and UMA. They require explicit opt-in and do not download checkpoints during ordinary test runs:

```bash
QUACC_RUN_OC25_INTEGRATION=1 QUACC_OC25_DEVICE=cuda \
    pytest tests/core/recipes/mlip_recipes/test_oc25.py
```

Set `QUACC_OC25_DEVICE=cpu` to use CPU inference. Optionally set `QUACC_OC25_STRUCTURE` to an ASE-readable structure file; otherwise the tests use a small Cu/water example intended only to exercise the integration. `QUACC_OC25_ESEN_CHECKPOINT` and `QUACC_OC25_UMA_CHECKPOINT` can select corresponding local checkpoint files. The tests establish integration parity and exercise a few optimization/MD steps; they do not establish DFT accuracy, relaxation convergence, or long-time dynamical stability.

With VASP and the recommended VASP 6.4 potentials configured, the separate
licensed-calculation smoke test can be run with:

```bash
QUACC_RUN_OC25_VASP=1 \
    pytest tests/core/recipes/vasp_recipes/jenkins/test_fairchem_oc25.py
```

It accepts the same optional `QUACC_OC25_STRUCTURE` file and otherwise uses a
small Cu/water interface. Ordinary test runs skip this calculation.
