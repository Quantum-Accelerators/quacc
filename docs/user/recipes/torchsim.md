# TorchSim Recipes

This tutorial introduces quacc's TorchSim recipes, which provide a high-level interface for running batched molecular dynamics simulations, geometry optimizations, and static calculations using machine learning potentials. The quacc TorchSim recipes wrap [TorchSim's high-level API](https://github.com/orbital-materials/torch-sim) while providing workflow engine compatibility, automatic trajectory handling, and standardized output schemas.

## Pre-Requisites

Before using the TorchSim recipes, you need to install TorchSim and at least one machine learning potential:

```bash
pip install torch-sim mace-torch
```

You should also be familiar with [ASE Atoms objects](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) and have read the [Intro to Jobs](recipes_intro.md) tutorial.

!!! Note

    Since we are not using a workflow engine for these examples, run the following in the command line:

    ```bash
    quacc set WORKFLOW_ENGINE None
    ```

## Overview

The quacc TorchSim recipes provide three main jobs:

1. **`relax_job`** - Geometry optimization using batched optimizers
2. **`md_job`** - Molecular dynamics simulations with various integrators
3. **`static_job`** - Single-point energy/force/stress calculations

Key features include:

- **Batch processing**: Simulate multiple systems in a single job
- **Autobatching**: Automatically manage GPU memory for large batches
- **Trajectory reporting**: Save trajectories and computed properties
- **Workflow compatibility**: Works with all quacc-supported workflow engines

## Basic Geometry Optimization

Let's start with a simple geometry optimization. We'll use a Lennard-Jones potential for demonstration, but the same approach works with any supported model type.

```python
from ase.build import bulk
import torch_sim as ts
from quacc.recipes.torchsim.core import relax_job
from quacc.schemas.torchsim import TSModelType

# Create an Argon FCC structure and perturb it
atoms = bulk("Ar", "fcc", a=5.26, cubic=True)
atoms.positions += 0.1  # Perturb to make optimization meaningful

# Run geometry optimization
result = relax_job(
    atoms=[atoms],
    model_type=TSModelType.LENNARD_JONES,
    model_path=None,
    optimizer=ts.Optimizer.fire,
    init_kwargs={"cell_filter": ts.CellFilter.unit},
    model_kwargs={"sigma": 3.405, "epsilon": 0.0104, "compute_stress": True},
)

# Access results
print(f"Final energy: {result['output']['energy'][0]:.6f} eV")
print(f"Optimized atoms: {result['atoms'][0]}")
```

### Understanding the Output

The `relax_job` returns a dictionary containing:

- **`atoms`**: List of optimized ASE Atoms objects
- **`output`**: Dictionary with `energy`, `forces`, and `stress` for each system
- **`dir_name`**: Path to the job directory containing trajectory files
- **`model_type`**, **`model_path`**: Model configuration used
- **`optimizer`**, **`convergence_fn`**: Optimization settings
- **`trajectory_reporter`**, **`autobatcher`**: Configuration details

## Molecular Dynamics

The `md_job` recipe runs molecular dynamics simulations. Here's an example using the Langevin thermostat for NVT dynamics:

```python
from ase.build import bulk
import torch_sim as ts
from quacc.recipes.torchsim.core import md_job
from quacc.schemas.torchsim import TSModelType

# Create a copper structure
atoms = bulk("Cu", "fcc", a=3.6, cubic=True)

# Run NVT molecular dynamics
result = md_job(
    atoms=[atoms],
    model_type=TSModelType.LENNARD_JONES,
    model_path=None,
    integrator=ts.Integrator.nvt_langevin,
    n_steps=100,
    temperature=300.0,  # Kelvin
    timestep=0.001,  # picoseconds
    model_kwargs={"sigma": 2.0, "epsilon": 0.1},
)

print(f"Final energy: {result['output']['energy'][0]:.6f} eV")
print(f"MD ran for {result['n_steps']} steps at {result['temperature']} K")
```

### Available Integrators

TorchSim provides several integrators accessible via `ts.Integrator`:

- `nvt_langevin` - Langevin thermostat (NVT)
- `nvt_nose_hoover` - Nosé-Hoover thermostat (NVT)
- `nve` - Microcanonical ensemble (NVE)
- `npt_langevin` - Langevin barostat (NPT)

## Static Calculations

For single-point calculations without system evolution, use `static_job`:

```python
from ase.build import bulk
from quacc.recipes.torchsim.core import static_job
from quacc.schemas.torchsim import TSModelType

# Create multiple structures
cu_atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
fe_atoms = bulk("Fe", "bcc", a=2.87, cubic=True)

# Run static calculations on both
result = static_job(
    atoms=[cu_atoms, fe_atoms],
    model_type=TSModelType.LENNARD_JONES,
    model_path=None,
    trajectory_reporter_dict={
        "prop_calculators": {1: ["potential_energy", "forces", "stress"]},
    },
    model_kwargs={"sigma": 2.5, "epsilon": 0.05, "compute_stress": True},
)

# Access energies for each system
for i, energy in enumerate(result["output"]["energy"]):
    print(f"System {i}: {energy:.6f} eV")
```

## Using Machine Learning Potentials

The TorchSim recipes support various machine learning potentials. Here's an example with MACE:

```python
from ase.build import bulk
import torch_sim as ts
from mace.calculators.foundations_models import download_mace_mp_checkpoint
from quacc.recipes.torchsim.core import relax_job
from quacc.schemas.torchsim import TSModelType

# Download and get path to MACE model
model_path = download_mace_mp_checkpoint("small")

# Create a silicon structure
atoms = bulk("Si", "diamond", a=5.43, cubic=True)
atoms.positions += 0.05  # Small perturbation

# Run optimization with MACE
result = relax_job(
    atoms=[atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    optimizer=ts.Optimizer.fire,
    convergence_fn="force",
    convergence_fn_kwargs={"force_tol": 0.01},  # eV/Å
    init_kwargs={"cell_filter": ts.CellFilter.unit},
)

print(f"Final energy: {result['output']['energy'][0]:.6f} eV")
```

### Supported Model Types

The `TSModelType` enum includes:

| Model Type | Description |
|------------|-------------|
| `MACE` | MACE models |
| `FAIRCHEM` | FairChem (EquiformerV2, eSCN) |
| `FAIRCHEMV1` | Legacy FairChem models |
| `MATTERSIM` | MatterSim models |
| `SEVENNET` | SevenNet models |
| `ORB` | Orb models |
| `GRAPHPESWRAPPER` | GraphPES models |
| `NEQUIPFRAMEWORK` | NequIP Framework models |
| `METATOMIC` | Metatomic models |
| `LENNARD_JONES` | Classical Lennard-Jones (testing) |

## Batch Processing Multiple Systems

One of TorchSim's key strengths is efficient batch processing. You can simulate multiple systems in parallel:

```python
from ase.build import bulk
import torch_sim as ts
from quacc.recipes.torchsim.core import relax_job
from quacc.schemas.torchsim import TSModelType

# Create multiple systems of different sizes
systems = [
    bulk("Cu", "fcc", a=3.6, cubic=True),
    bulk("Fe", "bcc", a=2.87, cubic=True),
    bulk("Cu", "fcc", a=3.6, cubic=True).repeat([2, 2, 2]),
    bulk("Fe", "bcc", a=2.87, cubic=True).repeat([2, 2, 2]),
]

# Perturb all structures
for atoms in systems:
    atoms.positions += 0.1

# Optimize all systems in a single job
result = relax_job(
    atoms=systems,
    model_type=TSModelType.LENNARD_JONES,
    model_path=None,
    optimizer=ts.Optimizer.fire,
    init_kwargs={"cell_filter": ts.CellFilter.unit},
    model_kwargs={"sigma": 2.5, "epsilon": 0.05, "compute_stress": True},
)

# Results for each system
for i, (atoms, energy) in enumerate(zip(result["atoms"], result["output"]["energy"])):
    print(f"System {i}: {len(atoms)} atoms, {energy:.6f} eV")
```

## Trajectory Reporting

Track simulation progress and save trajectories using `trajectory_reporter_dict`:

```python
from ase.build import bulk
import torch_sim as ts
from quacc.recipes.torchsim.core import md_job
from quacc.schemas.torchsim import TSModelType

atoms = bulk("Cu", "fcc", a=3.6, cubic=True)

# Configure trajectory reporting
trajectory_reporter = {
    "state_frequency": 10,  # Save state every 10 steps
    "prop_calculators": {
        5: ["potential_energy", "kinetic_energy", "temperature"],
    },
}

result = md_job(
    atoms=[atoms],
    model_type=TSModelType.LENNARD_JONES,
    model_path=None,
    integrator=ts.Integrator.nvt_langevin,
    n_steps=100,
    temperature=500.0,
    timestep=0.001,
    trajectory_reporter_dict=trajectory_reporter,
    model_kwargs={"sigma": 2.0, "epsilon": 0.1},
)

# Trajectory files are saved in result["dir_name"]
print(f"Trajectory saved to: {result['dir_name']}")
print(f"Reporter config: {result['trajectory_reporter']}")
```

### Available Property Calculators

Properties that can be tracked during simulations:

- `potential_energy` - System potential energy
- `forces` - Atomic forces
- `stress` - System stress tensor
- `kinetic_energy` - Kinetic energy (MD only)
- `temperature` - Instantaneous temperature (MD only)

## Autobatching

For large batches that exceed GPU memory, enable autobatching:

```python
from ase.build import bulk
import torch_sim as ts
from quacc.recipes.torchsim.core import relax_job
from quacc.schemas.torchsim import TSModelType

# Create many systems
systems = [bulk("Cu", "fcc", a=3.6, cubic=True).repeat([2, 2, 2]) for _ in range(20)]
for atoms in systems:
    atoms.positions += 0.1

# Enable autobatching with custom settings
autobatcher = {
    "memory_scales_with": "n_atoms",
    "max_memory_scaler": 260,  # Adjust based on GPU memory
}

result = relax_job(
    atoms=systems,
    model_type=TSModelType.LENNARD_JONES,
    model_path=None,
    optimizer=ts.Optimizer.fire,
    autobatcher_dict=autobatcher,
    init_kwargs={"cell_filter": ts.CellFilter.unit},
    model_kwargs={"sigma": 2.5, "epsilon": 0.05, "compute_stress": True},
)

print(f"Optimized {len(result['atoms'])} systems")
```

You can also set `autobatcher_dict=True` for automatic configuration.

## Convergence Criteria

For geometry optimization, two convergence criteria are available:

### Force Convergence (Default)

```python
result = relax_job(
    atoms=[atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    optimizer=ts.Optimizer.fire,
    convergence_fn="force",
    convergence_fn_kwargs={"force_tol": 0.01},  # eV/Å
    init_kwargs={"cell_filter": ts.CellFilter.unit},
)
```

### Energy Convergence

```python
result = relax_job(
    atoms=[atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    optimizer=ts.Optimizer.fire,
    convergence_fn="energy",
    convergence_fn_kwargs={"energy_tol": 1e-6},  # eV
)
```

## Cell Optimization

To optimize both atomic positions and cell parameters, use the `cell_filter` option:

```python
from ase.build import bulk
import torch_sim as ts
from quacc.recipes.torchsim.core import relax_job
from quacc.schemas.torchsim import TSModelType

atoms = bulk("Si", "diamond", a=5.5)  # Slightly wrong lattice constant
atoms.positions += 0.05

result = relax_job(
    atoms=[atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    optimizer=ts.Optimizer.fire,
    init_kwargs={"cell_filter": ts.CellFilter.unit},  # Enable cell optimization
)
```

Available cell filters via `ts.CellFilter`:

- `unit` - Full cell relaxation
- `scalar` - Isotropic volume relaxation
- `off` - Fixed cell (positions only)

## Complete Workflow Example

Here's a complete example combining multiple jobs:

```python
from ase.build import bulk
import torch_sim as ts
from mace.calculators.foundations_models import download_mace_mp_checkpoint
from quacc.recipes.torchsim.core import relax_job, md_job, static_job
from quacc.schemas.torchsim import TSModelType

# Setup
model_path = download_mace_mp_checkpoint("small")
atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat([2, 2, 2])
atoms.positions += 0.1

# Step 1: Optimize the structure
relax_result = relax_job(
    atoms=[atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    optimizer=ts.Optimizer.fire,
    init_kwargs={"cell_filter": ts.CellFilter.unit},
)
optimized_atoms = relax_result["atoms"][0]
print(f"Relaxation energy: {relax_result['output']['energy'][0]:.4f} eV")

# Step 2: Run MD to equilibrate
md_result = md_job(
    atoms=[optimized_atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    integrator=ts.Integrator.nvt_langevin,
    n_steps=500,
    temperature=300.0,
    timestep=0.002,
    trajectory_reporter_dict={
        "state_frequency": 50,
        "prop_calculators": {10: ["potential_energy", "temperature"]},
    },
)
equilibrated_atoms = md_result["atoms"][0]
print(f"MD final energy: {md_result['output']['energy'][0]:.4f} eV")

# Step 3: Static calculation for final properties
static_result = static_job(
    atoms=[equilibrated_atoms],
    model_type=TSModelType.MACE,
    model_path=model_path,
    trajectory_reporter_dict={
        "prop_calculators": {1: ["potential_energy", "forces", "stress"]},
    },
)
print(f"Final energy: {static_result['output']['energy'][0]:.4f} eV")
print(f"Final forces shape: {len(static_result['output']['forces'][0])} atoms")
```

## Comparison with TorchSim Direct API

The quacc recipes wrap TorchSim's high-level functions with additional features:

| Feature | TorchSim Direct | quacc Recipe |
|---------|-----------------|--------------|
| Batch processing | `ts.optimize(system=systems, ...)` | `relax_job(atoms=systems, ...)` |
| Workflow engines | Manual integration | Built-in support |
| Result schemas | Custom handling | Standardized dictionaries |
| Trajectory paths | Manual management | Automatic job directory |
| Model loading | Manual instantiation | `model_type` + `model_path` |

## Conclusion

The quacc TorchSim recipes provide a powerful interface for running batched atomistic simulations with machine learning potentials. Key advantages include:

1. **Workflow compatibility** - Works with Prefect, Parsl, Dask, and other engines
2. **Batch efficiency** - Process multiple systems in parallel on GPU
3. **Automatic memory management** - Autobatching handles large datasets
4. **Standardized outputs** - Consistent result schemas across all jobs
5. **Trajectory management** - Automatic file handling in job directories

For more details on the underlying TorchSim functionality, see the [TorchSim documentation](https://orbital-materials.github.io/torch-sim/).
