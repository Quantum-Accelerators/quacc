# Psi4

🚧 Under Construction 🚧

Recipes based on the [Psi4](https://psicode.org/) code.

## Available Recipes

| Recipe                                                                                                                                           | Description                       | Calculator Parameters                                                                                                                                                                 | Other Parameters |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [`psi4.core.static_job`](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/psi4/core.html#quacc.recipes.psi4.core.static_job) | Static calculation on a molecule. | `#!Python {"mem": "16GB", "num_threads": "max", "method": method, "basis": basis, "charge": charge, "multiplicity": multiplicity, "reference": "uks" if multiplicity > 1 else "rks"}` | N/A              |

## Examples

### Core Recipes

```mermaid
graph LR
  A[Input] --> B(Psi4 Static) --> C[Output];
```

=== "No Workflow Engine"

    ```python
    from ase.build import molecule
    from quacc.recipes.psi4.core import static_job

    atoms = molecule("O2")
    result = static_job(atoms, charge=0, multiplicity=3, method="wb97m-v", basis="def2-svp")
    ```

=== "Covalent"

    ```python
    import covalent as ct
    from ase.build import molecule
    from quacc import flow
    from quacc.recipes.psi4.core import static_job

    workflow = flow(static_job)
    atoms = molecule("O2")

    dispatch_id = workflow(
        atoms, charge=0, multiplicity=3, method="wb97m-v", basis="def2-svp"
    )
    result = ct.get_result(dispatch_id, wait=True)
    ```

=== "Parsl"

    ```python
    from ase.build import molecule
    from quacc.recipes.psi4.core import static_job

    atoms = molecule("O2")
    future = static_job(atoms, charge=0, multiplicity=3, method="wb97m-v", basis="def2-svp")
    result = future.result()
    ```
