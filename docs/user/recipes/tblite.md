# TBLite

Recipes based on the [TBLite](https://github.com/tblite/tblite) code, which has a Python interface to the [xtb](https://github.com/grimme-lab/xtb) library for semi-empirical tight-binding calculations.

## Available Recipes

| Recipe                                                                                                                                                 | Description                                     | Calculator Parameters | Other Parameters                                                |
| ------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------- | --------------------- | --------------------------------------------------------------- |
| [`tblie.core.relax_job`](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/tblite/core.html#quacc.recipes.tblite.core.relax_job)    | Relax a molecule or crystal.                    | `#!Python {}`         | `#!Python {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}` |
| [`tblite.core.static_job`](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/tblite/core.html#quacc.recipes.tblite.core.static_job) | Static calculation on a molecule or crystal.    | `#!Python {}`         | N/A                                                             |
| [`tblite.core.freq_job`](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/tblite/core.html#quacc.recipes.tblite.core.freq_job)     | Frequency calculation on a molecule or crystal. | `#!Python {}`         | N/A                                                             |

## Examples

### Core Recipes

```mermaid
graph LR
  A[Input] --> B(GFN2-xTB Cell Relax) --> C[Output];
```

=== "No Workflow Engine"

    ```python
    from ase.build import bulk
    from quacc.recipes.tblite.core import relax_job

    atoms = bulk("C")
    result = relax_job(atoms, relax_cell=True)
    ```

=== "Covalent"

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.tblite.core import relax_job

    atoms = bulk("C")
    dispatch_id = flow(relax_job)(atoms, relax_cell=True)
    result = ct.get_result(dispatch_id, wait=True)
    ```

=== "Parsl"

    ```python
    from ase.build import bulk
    from quacc.recipes.tblite.core import relax_job

    atoms = bulk("C")
    future = relax_job(atoms, relax_cell=True)
    result = future.result()
    ```

```mermaid
graph LR
  A[Input] --> B(GFN1-xTB Static) --> C[Output];
```

=== "No Workflow Engine"

    ```python
    from ase.build import bulk
    from quacc.recipes.tblite.core import static_job

    atoms = bulk("C")
    result = static_job(atoms, method="GFN1-xTB")
    ```

=== "Covalent"

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.tblite.core import static_job

    atoms = bulk("C")
    dispatch_id = flow(static_job)(atoms, method="GFN1-xTB")
    result = ct.get_result(dispatch_id, wait=True)
    ```

=== "Parsl"

    ```python
    from ase.build import bulk
    from quacc.recipes.tblite.core import static_job

    atoms = bulk("C")
    future = static_job(atoms, method="GFN1-xTB")
    result = future.result()
    ```

```mermaid
graph LR
  A[Input] --> B(GFN2-xTB Frequency Analysis) --> C[Output];
```

=== "No Workflow Engine"

    ```python
    from ase.build import molecule
    from quacc.recipes.tblite.core import freq_job

    atoms = molecule("N2")
    result = freq_job(atoms)
    ```

=== "Covalent"

    ```python
    import covalent as ct
    from ase.build import molecule
    from quacc import flow
    from quacc.recipes.tblite.core import freq_job

    atoms = molecule("N2")
    dispatch_id = flow(freq_job)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    ```

=== "Parsl"

    ```python
    from ase.build import molecule
    from quacc.recipes.tblite.core import freq_job

    atoms = molecule("N2")
    future = freq_job(atoms)
    result = future.result()
    ```
